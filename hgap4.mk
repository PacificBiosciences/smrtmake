# SGE queue name to submit jobs to
QUEUE ?= huasm
# Size of the genome
GENOME_SIZE ?= 5000000
# Splits data into this many chunks, each chunk processed independently
CHUNK_SIZE ?= 3
# How many threads a process will use (also how many SGE slots will be requested)
NPROC ?= 15
# Local temp root directory, must have write access and a decent amount of space (~100GB)
LOCALTMP ?= /scratch

SPLITBESTN = $(shell echo $$(( 30/$(CHUNK_SIZE)+2 )))

SHELL = /bin/bash
TASKDIRS = filter correct assemble polish log

SMRTETC = $(SEYMOUR_HOME)/analysis/etc
PWD = $(shell pwd)

CHUNKS := $(shell seq 1 $(CHUNK_SIZE))

# Chunks are maintained throughout workflow. Need to add some elasticity to handle
# larger datasets
BAXFOFNS := $(foreach c,$(CHUNKS),$(shell printf "input.chunk%03dof%03d.fofn" $(c) $(CHUNK_SIZE)))
REGFOFNS := $(BAXFOFNS:input.%=filter/regions.%)
SUBLENGTHS := $(BAXFOFNS:input.%.fofn=filter/subreads.%.lengths)
LONGFASTA := $(BAXFOFNS:input.%.fofn=filter/longreads.%.fasta)
MAPPEDM4 := $(BAXFOFNS:input.%.fofn=correct/seeds.%.m4)
MAPPEDM4FIX := $(BAXFOFNS:input.%.fofn=correct/seeds.%.m4.fix)
MAPPEDM4FILT := $(BAXFOFNS:input.%.fofn=correct/seeds.%.m4.filt)
CORRECTED := $(BAXFOFNS:input.%.fofn=correct/corrected.%.fasta)
CMPH5 := $(BAXFOFNS:input.%.fofn=polish/aligned_reads.%.cmp.h5)

CUTOFF = filter/longreads.cutoff
ALLREGFOFN = filter/regions.fofn
ALLBAXFOFN = input.fofn

QSUB = qsub -S /bin/bash -cwd -b y -sync y -V -q $(QUEUE) -e $(PWD)/log/ -o $(PWD)/log/

all : assembly

prepare : 
	mkdir -p $(TASKDIRS)

assembly : polish/polished_assembly.fasta | prepare

## Assembly polishing ##

polish/polished_assembly.fasta : polish/aligned_reads.cmp.h5 polish/reference cmph5_chemistry_loaded
	$(QSUB) -N polish -pe smp $(NPROC) variantCaller.py -P $(SMRTETC)/algorithm_parameters/2014-03 \
	-v -j $(NPROC) --algorithm=quiver $< -r $(word 2,$^)/sequence/reference.fasta -o polish/corrections.gff \
	-o $@ -o $(@:.fasta=.fastq.gz)

cmph5_chemistry_loaded : filter/chemistry_mapping.xml polish/aligned_reads.cmp.h5
	loadSequencingChemistryIntoCmpH5.py --xml $< --h5 $(word 2,$^) && touch chemistry_loaded

polish/aligned_reads.cmp.h5 : $(CMPH5)
	assertCmpH5NonEmpty.py --debug $^
	$(QSUB) -N mergesort 'cmph5tools.py -vv merge --outFile=$@ $^ && cmph5tools.py -vv sort --deep --inPlace $@'
	#h5repack -f GZIP=1 $@ $@_TMP && mv $@_TMP $@

cmph5 : $(CMPH5) ;

$(CMPH5) : polish/aligned_reads.%.cmp.h5 : input.%.fofn filter/regions.%.fofn | polish/reference
	$(QSUB) -N res.$* -pe smp $(NPROC) "pbalign.py $< $| $@ --seed=1 --minAccuracy=0.75 --minLength=50 --algorithmOptions=\"-useQuality -minMatch 12 -bestn 10 -minPctIdentity 70.0\" --hitPolicy=randombest --tmpDir=$(LOCALTMP) -vv --nproc=$(NPROC) --regionTable=$(word 2,$^) && loadPulses $< $@ -metrics DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag -byread"

polish/reference : assemble/draft_assembly.fasta
	referenceUploader --skipIndexUpdate -c -n "reference" -p polish -f $<  --saw="sawriter -blt 8 -welter" --samIdx="samtools faidx"

##

## Assembly draft ##
draft : assemble/draft_assembly.fasta ;

assemble/draft_assembly.fasta : assemble/utg.finished
	tigStore -g assemble/celera-assembler.gkpStore -t assemble/celera-assembler.tigStore 1 -d properties -U \
	| awk 'BEGIN{t=0}$$1=="numFrags"{if($$2>1){print t, $$2}t++}' | sort -nrk2,2 > assemble/unitig.lst
	$(QSUB) -N draft -pe smp $(NPROC) 'tmp=$$(mktemp -d -p $(LOCALTMP)); tmp=$$tmp cap=$(PWD)/assemble/celera-assembler utg=$(PWD)/assemble/unitig.lst cns=$(PWD)/$@ nproc=$(NPROC) pbutgcns_wf.sh'

assemble/utg.finished : assemble/utg.spec assemble/utg.frg
	$(QSUB) -pe smp $(NPROC) runCA -d assemble -p celera-assembler -s $^
	touch $@

assemble/utg.spec : correct/corrected.fasta
	runCASpecWriter.py  -vv --bitTable=$(SMRTETC)/celeraAssembler/bitTable \
	--interactiveTmpl=$(SMRTETC)/cluster/SGE/interactive.tmpl \
	--smrtpipeRc=$(SMRTETC)/smrtpipe.rc --genomeSize=$(GENOME_SIZE) --defaultFrgMinLen=500 \
	--xCoverage=20 --ovlErrorRate=0.06 --ovlMinLen=40 --merSize=14 --corrReadsFasta=$< \
	--specOut=$@ --sgeName=utg --gridParams="useGrid:0, scriptOnGrid:0, frgCorrOnGrid:0, ovlCorrOnGrid:0" \
	--maxSlotPerc=1 $(SMRTETC)/celeraAssembler/unitig.spec

assemble/utg.frg : $(CORRECTED)
	fastqToCA -technology sanger -type sanger -libraryname corr $(patsubst %,-reads %,$(^:.fasta=.fastq)) > $@

correct/corrected.fasta : $(CORRECTED)
	cat $^ > $@

##

## Correction (optimizations available here) ##
correction : $(CORRECTED) ;

$(CORRECTED) : correct/corrected.%.fasta : correct/seeds.%.m4.filt correct/seeds.m4.fofn
	$(QSUB) -N corr.$* -pe smp $(NPROC) 'tmp=$$(mktemp -d -p $(LOCALTMP)); mym4=$(PWD)/$< allm4=$(PWD)/$(word 2,$^) basfofn=$(PWD)/$(ALLBAXFOFN) rgnfofn=$(PWD)/$(ALLREGFOFN) bestn=24 nproc=$(NPROC) fasta=$(PWD)/$@ fastq=$(PWD)/$(@:.fasta=.fastq) tmp=$$tmp pbdagcon_wf.sh; rm -rf $$tmp'

correct/seeds.m4.fofn : $(MAPPEDM4FILT)
	echo $(^:%=$(PWD)/%) | sed 's/ /\n/g' > $@

$(MAPPEDM4FILT) : correct/seeds.%.m4.filt : correct/seeds.%.m4.fix
	filterm4.py $< > $@

##

## Read overlap using BLASR ##
overlap : $(MAPPEDM4FIX)

$(MAPPEDM4FIX) : correct/seeds.%.m4.fix : correct/seeds.%.m4
	perl -F'\s' -ane '@c=$$F[0]=~m#(\d+)_(\d+)$$#;$$F[5]-=$$c[0];$$F[6]-=$$c[0];$$F[7]=$$c[1]-$$c[0];print join " ",@F;print "\n"' $< > $@ 

$(MAPPEDM4) : correct/seeds.%.m4 : filter/longreads.%.fasta $(ALLREGFOFN)
	$(QSUB) -N blasr.$* -pe smp $(NPROC) blasr input.fofn $< -out $@ -regionTable $(word 2,$^) -m 4 -nproc $(NPROC) -bestn $(SPLITBESTN) -nCandidates $(SPLITBESTN) -maxScore -1000 -maxLCPLength 16 -minMatch 14

$(ALLREGFOFN) : $(REGFOFNS)
	sed 's;.*Results/\([^\.]*\.[1-3]\).*;filter/\1.rgn.h5;' input.fofn > $@
##

## Generating the long seed reads for mapping ##
longreads : $(LONGFASTA);

# bash5minlen.py can be found in utils/
$(LONGFASTA) : filter/longreads.%.fasta : input.%.fofn filter/regions.%.fofn $(CUTOFF)
	bash5minlen.py $^ > $@

$(CUTOFF) : $(SUBLENGTHS)
	sort -nrmk1,1 $^ | awk '{t+=$$1;if(t>=$(GENOME_SIZE)*30){print $$1;exit;}}' > $@
	
# bash5lengths.py can be found in utils/
$(SUBLENGTHS) : filter/subreads.%.lengths : input.%.fofn filter/regions.%.fofn
	bash5lengths.py $^ > $@
##

## Filtering ##
regions : $(REGFOFNS) ;

$(REGFOFNS) : filter/regions.%.fofn : input.%.fofn | prepare
	$(QSUB) -N filt.$* filter_plsh5.py --filter='MinReadScore=0.80,MinSRL=500,MinRL=100' \
	--trim='True' --outputDir=filter --outputFofn=$@ $<

filter/chemistry_mapping.xml : filter/movie_metadata
	makeChemistryMapping.py --metadata=filter/movie_metadata --output=filter/chemistry_mapping.xml --mapping_xml=$(SMRTETC)/algorithm_parameters/2013-09/mapping.xml

filter/movie_metadata : input.fofn
	mkdir -p $@
	sed 's#\(.*\)Analysis_Results/\(m[^\.]*\).*#\1\2.metadata.xml#' $< | sort -u | xargs ln -s -t $@

##

## Initial chunking ##
inputs : $(BAXFOFNS) ;

$(BAXFOFNS) : input.fofn
	awk 'BEGIN{c=1}{print $$0 > sprintf("input.chunk%03dof%03d.fofn", c++%$(CHUNK_SIZE)+1, $(CHUNK_SIZE))}' $<
##

clean :
	rm -rf $(TASKDIRS)
	rm -f chemistry_loaded
	rm -f input.chunk*
