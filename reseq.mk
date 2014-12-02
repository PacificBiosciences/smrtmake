include chunk.mk

NPROC ?= 8
LOCALTMP ?= /scratch
QUEUE ?= secondary
QSUB = qsub -S /bin/bash -cwd -b y -sync y -V -q $(QUEUE) -j y -e $(PWD)/log/ -o $(PWD)/log/
REFERENCE ?= reference.fasta
REFERENCE_IDX := $(REFERENCE).fai

REGFOFN_CHUNKS := $(BAXFOFN_CHUNKS:input.%.fofn=filter/regions.%.fofn)
CMPH5_CHUNKS := $(BAXFOFN_CHUNKS:input.%.fofn=mapping/aligned_reads.%.cmp.h5)

TASKDIRS = filter mapping quiver log

ALIGN = pbalign $< $| $@ --seed=1 --minAccuracy=0.75 --minLength=50 \
	--algorithmOptions=\"-useQuality -minMatch 12 -bestn 10 -minPctIdentity 70.0\" \
	--hitPolicy=randombest --tmpDir=$(LOCALTMP) -vv --nproc=$(NPROC) --regionTable=$(word 2,$^) \
	&& loadPulses $< $@ -metrics DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag -byread

QUIVER = variantCaller.py -v -j $(NPROC) --algorithm=quiver $< -r $(word 2,$^)\
	 -o quiver/corrections.gff -o $@ -o $(@:.fasta=.fastq.gz)

all : quiver

prepare : 
	mkdir -p $(TASKDIRS)

quiver : quiver/results.fasta ;

quiver/results.fasta : mapping/aligned_reads.cmp.h5 $(REFERENCE) $(REFERENCE_IDX) cmph5_chemistry_loaded
	$(QSUB) -N polish -pe smp $(NPROC) $(QUIVER)

cmph5_chemistry_loaded : $(BAXFOFN) mapping/aligned_reads.cmp.h5
	loadChemistry.py  $^ && touch $@

mapping : mapping/aligned_reads.cmp.h5 ;

mapping/aligned_reads.cmp.h5 : $(CMPH5_CHUNKS)
	assertCmpH5NonEmpty.py --debug $^
	$(QSUB) -N mergesort 'cmph5tools.py -vv merge --outFile=$@ $^;cmph5tools.py -vv sort --deep --inPlace $@'
	h5repack -f GZIP=1 $@ $@_TMP && mv $@_TMP $@

$(CMPH5_CHUNKS) : mapping/aligned_reads.%.cmp.h5 : input.%.fofn filter/regions.%.fofn | $(REFERENCE)
	$(QSUB) -N res.$* -pe smp $(NPROC) "$(ALIGN)"

$(REGFOFN_CHUNKS) : filter/regions.%.fofn : input.%.fofn | prepare
	$(QSUB) -N filt.$* filter_plsh5.py --filter='MinReadScore=0.80,MinSRL=500,MinRL=100' \
	--trim='True' --outputDir=filter --outputFofn=$@ $<

$(REFERENCE_IDX) : $(REFERENCE) 
	samtools faidx $<

$(REFERENCE) : ;

clean : 
	rm -rf $(TASKDIRS)
	rm -f cmph5_chemistry_loaded
	rm -f $(BAXFOFN_CHUNKS) 
	rm -f $(INIT_CHUNKS) 
