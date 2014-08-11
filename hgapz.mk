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

INPUT ?= input.fofn

SHELL = /bin/bash
TASKDIRS = filter correct assemble polish log
SUBPFX = subreads
DALPFX = $(SUBPFX).$(SUBPFX)

SMRTETC = $(SEYMOUR_HOME)/analysis/etc
PWD = $(shell pwd)

VPATH := $(shell sed 's:^\(.*\)/.*:\1:' $(INPUT) | sort -u)
MOVIES := $(shell sed 's:.*/\([^.]*\).[1-3].bax.h5:\1:' input.fofn | sort -u)
SUBFASTA := $(MOVIES:%=filter/%.$(SUBPFX).fasta)
SUBLENGTHS := $(SUBFASTA:%.fasta=%.lengths)
SUBSEEDS := $(SUBLENGTHS:%.lengths=%.seeds)
DALFILES = $(DALPFX).C0.las $(DALPFX).C1.las $(DALPFX).C2.las $(DALPFX).C3.las \
           $(DALPFX).N0.las $(DALPFX).N1.las $(DALPFX).N2.las $(DALPFX).N3.las

DAZALN := $(DALFILES:%=correct/%)
DAZALNSORT := $(DAZALN:%.las=%.S.las)

# Chunks are maintained throughout workflow. Need to add some elasticity to handle
# larger datasets
CHUNKS := $(shell seq 1 $(CHUNK_SIZE))
SPLITBESTN = $(shell echo $$(( 30/$(CHUNK_SIZE)+2 )))
BAXFOFNS := $(foreach c,$(CHUNKS),$(shell printf "input.chunk%03dof%03d.fofn" $(c) $(CHUNK_SIZE)))
REGFOFNS := $(BAXFOFNS:input.%=filter/regions.%)
DAZ2M4 := $(BAXFOFNS:input.%.fofn=correct/daz2m4.%.seeds)
MAPPEDM4 := $(BAXFOFNS:input.%.fofn=correct/seeds.%.m4)
MAPPEDM4FILT := $(BAXFOFNS:input.%.fofn=correct/seeds.%.m4.filt)
CORRECTED := $(BAXFOFNS:input.%.fofn=correct/corrected.%.fasta)
CMPH5 := $(BAXFOFNS:input.%.fofn=polish/aligned_reads.%.cmp.h5)

CUTOFF = filter/longreads.cutoff
QUERYFOFN = filter/subreads.fofn

QSUB = qsub -S /bin/bash -cwd -b y -sync y -V -q $(QUEUE) -e $(PWD)/log/ -o $(PWD)/log/

all : assembly

prepare : 
	mkdir -p $(TASKDIRS)

assembly : polish/polished_assembly.fasta | prepare

## Assembly polishing ##

polish/polished_assembly.fasta : polish/aligned_reads.cmp.h5 polish/reference 
	$(QSUB) -N polish -pe smp $(NPROC) variantCaller.py -P $(SMRTETC)/algorithm_parameters/2014-03 \
	-v -j $(NPROC) --algorithm=quiver $< -r $(word 2,$^)/sequence/reference.fasta -o polish/corrections.gff \
	-o $@ -o $(@:.fasta=.fastq.gz)

polish/aligned_reads.cmp.h5 : $(CMPH5)
	assertCmpH5NonEmpty.py --debug $^
	$(QSUB) -N mergesort 'cmph5tools.py -vv merge --outFile=$@ $^;cmph5tools.py -vv sort --deep --inPlace $@'
	#h5repack -f GZIP=1 $@ $@_TMP && mv $@_TMP $@

$(CMPH5) : polish/aligned_reads.%.cmp.h5 : input.%.fofn filter/regions.%.fofn | polish/reference
	$(QSUB) -N res.$* -pe smp $(NPROC) "pbalign $< $| $@ --forQuiver --seed=1 --minAccuracy=0.75 --minLength=50 --algorithmOptions=\"-useQuality -minMatch 12 -bestn 10 -minPctIdentity 70.0\" --hitPolicy=randombest --tmpDir=$(LOCALTMP) -vv --nproc=$(NPROC) --regionTable=$(word 2,$^) && loadPulses $< $@ -metrics DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag -byread"

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
	runCASpecWriter.py  --bitTable=$(SMRTETC)/celeraAssembler/bitTable \
	--interactiveTmpl=$(SMRTETC)/cluster/SGE/interactive.tmpl \
	--smrtpipeRc=$(SMRTETC)/smrtpipe.rc --genomeSize=$(GENOME_SIZE) --defaultFrgMinLen=500 \
	--xCoverage=20 --ovlErrorRate=0.06 --ovlMinLen=40 --merSize=14 --corrReadsFasta=$< \
	--specOut=$@ --sgeName=utg --gridParams="useGrid:1, scriptOnGrid:1, frgCorrOnGrid:1, ovlCorrOnGrid:1" \
	--maxSlotPerc=1 $(SMRTETC)/celeraAssembler/unitig.spec

assemble/utg.frg : $(CORRECTED)
	fastqToCA -technology sanger -type sanger -libraryname corr $(patsubst %,-reads %,$(^:.fasta=.fastq)) > $@

correct/corrected.fasta : $(CORRECTED)
	cat $^ > $@

##

## Correction ##
correct : $(CORRECTED) ;

$(CORRECTED) : correct/corrected.%.fasta : correct/seeds.%.m4.filt correct/seeds.m4.fofn filter/subreads.fasta
	$(QSUB) -N corr.$* -pe smp $(NPROC) 'tmp=$$(mktemp -d -p $(LOCALTMP)); mym4=$(PWD)/$< allm4=$(PWD)/$(word 2,$^) subreads=$(PWD)/$(word 3, $^) bestn=24 cov=6 nproc=$(NPROC) fasta=$(PWD)/$@ fastq=$(PWD)/$(@:.fasta=.fastq) tmp=$$tmp pbdagcon_wf.sh; rm -rf $$tmp'

correct/seeds.m4.fofn : $(MAPPEDM4FILT)
	echo $(^:%=$(PWD)/%) | sed 's/ /\n/g' > $@

filter : $(MAPPEDM4FILT) ;

$(MAPPEDM4FILT) : correct/seeds.%.m4.filt : correct/seeds.%.m4 
	filterm4.py $< > $@

# Minor edit to header and uppercase the bases
filter/subreads.fasta : $(SUBFASTA)
	sed '/^[^>]/ y/acgt/ACGT/;s/ RQ=.*//' $^ > $@
##

## The key bridgepoint: adapt daligner output to m4 format ##
adapt : $(MAPPEDM4) ;

# Yeah, I wrote it in perl ... so what?
$(MAPPEDM4) : correct/seeds.%.m4 : correct/daz2m4.%.seeds correct/$(DALPFX).merge.las $(SUBLENGTHS)
	daz2m4.pl $^ > $@

$(DAZ2M4) : rechunk.done ;

# Arrange the alignment data for processing into HGAP chunks. 
# daz2m4.<chunk>.seeds
rechunk.done : $(SUBSEEDS)
	awk '{print > sprintf("correct/daz2m4.chunk%03dof%03d.seeds", ++c%$(CHUNK_SIZE)+1, $(CHUNK_SIZE))}' $^
	@touch rechunk.done

##

## Read overlap using DALIGNER
# XXX: Need to see how blocks affect this output
dazaln : correct/$(DALPFX).merge.las ;

correct/$(DALPFX).merge.las : $(DAZALNSORT)
	LAmerge $@ $^

$(DAZALNSORT) : correct/%.S.las : correct/%.las correct/daligner.done
	LAsort $<

correct/daligner.done : correct/$(SUBPFX).db
	$(QSUB) -N daligner -pe smp 5 daligner -t9 $< $<
	mv $(DALFILES) correct/
	touch $@

$(DAZALN) : ;

##

## Find seed read cutoff and generate seeds
seeds : $(SUBSEEDS) ;

$(SUBSEEDS) : %.seeds : $(CUTOFF) %.lengths 
	awk -v cut=$$(cat $<) '($$1 >= cut){print}' $(lastword $^) > $@

$(CUTOFF) : $(SUBLENGTHS)
	sort -nrmk1,1 $^ | awk '{t+=$$1;if(t>=$(GENOME_SIZE)*30){print $$1;exit;}}' > $@

# NOTE: we this also emits what should be daz id mapping table
$(SUBLENGTHS) : filter/%.$(SUBPFX).lengths : filter/%.$(SUBPFX).fasta
	fastalength $< | awk '{print $$0, ++c}' | sort -nrk1,1 > $@
##

## Integration with dazzler flow, this has some explicit block enforcements
#  limited by the size of the daz internal read id.  All movies must be loaded
#  independently.  After creating the database, it must be split into blocks,
#  or chunks, if the # of reads is greater than uint16 (65,536).  This block
#  scheme will not be the same as the HGAP chunk scheme, so we'll have to 
#  manage dazzler blocks with HGAP chunks.
# XXX: Test on large enough dataset to necessitate block splitting.
dazdb : correct/$(SUBPFX).db

# NOTE: the order here determines how the daz internal ids get asseigned in the 
# database, important bridge for jumping between the HGAP/Dazzler worlds.
correct/$(SUBPFX).db : $(SUBFASTA)
	fasta2DB $@ $^

subreads : $(SUBFASTA) ;

$(SUBFASTA) : filter/%.$(SUBPFX).fasta : %.1.bax.h5 %.2.bax.h5 %.3.bax.h5 | prepare
	dextract -s800 $^ > $@

##

## Standard region filtering, still needed for quiver ##
$(REGFOFNS) : filter/regions.%.fofn : input.%.fofn | prepare
	$(QSUB) -N filt.$* filter_plsh5.py --filter='MinReadScore=0.80,MinSRL=500,MinRL=100' \
	--trim='True' --outputDir=filter --outputFofn=$@ $<

##

## Initial chunking (for HGAP-related workflow steps) ##
$(BAXFOFNS) : chunkinput.done ;

chunkinput.done : input.fofn
	awk 'BEGIN{c=1}{print $$0 > sprintf("input.chunk%03dof%03d.fofn", c++%$(CHUNK_SIZE)+1, $(CHUNK_SIZE))}' $<
	@touch $@
##

clean :
	rm -rf $(TASKDIRS)
	rm -f input.chunk*
	rm -f $(DALPFX)*
	rm -f rechunk.done
	rm -f chunkinput.done
