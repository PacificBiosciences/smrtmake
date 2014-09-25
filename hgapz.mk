# SGE queue name to submit jobs to
QUEUE ?= huasm
# Size of the genome
GENOME_SIZE ?= 3120000000
# Splits data into this many chunks, each chunk processed independently
CHUNK_SIZE ?= 3
# How many threads a process will use (also how many SGE slots will be requested)
NPROC ?= 15
# Local temp root directory, must have write access and a decent amount of space (~100GB)
LOCALTMP ?= /scratch

SHELL = /bin/bash

INPUT ?= input.fofn
SUBPFX ?= s
BLOCKD = blks
TASKDIRS = filter correct assemble polish log $(BLOCKD)

SMRTETC = $(SEYMOUR_HOME)/analysis/etc
PWD = $(shell pwd)

VPATH := $(shell sed 's:^\(.*\)/.*:\1:' $(INPUT) | sort -u)
MOVIES := $(shell sed 's:.*/\([^.]*\).[1-3].bax.h5:\1:' $(INPUT) | sort -u)
BAXFILES := $(shell cat $(INPUT))
SUBFASTA := $(MOVIES:%=filter/%.$(SUBPFX).fasta)
SUBLENGTHS := $(SUBFASTA:%.fasta=%.lengths)

# Dazzler block logic.  Need to know the number of blocks in order to construct
# this properly.  Creates 2TN^2 .las files, where T = 4 threads and N = number 
# of blocks. For perspective, arabidopsis will create 14,112 .las files.

DABLOCKS ?= 1

# Daligner output format: <dbname>.<block a>.<dbname>.<block b>.<strand><thread>

BLOCKFMT = $(SUBPFX).%d
BLKDONEFMT = $(BLOCKD)/$(BLOCKFMT)/$(BLOCKFMT)_$(BLOCKFMT).done
#DALDONE := $(shell python -c 'print " ".join(["$(BLKDONEFMT)" % (i, i, j) for i in xrange(1,$(DABLOCKS)+1) for j in xrange(i, $(DABLOCKS)+1)])')
DALDONE := $(shell python -c 'print " ".join(["$(BLKDONEFMT)" % ($(BLOCK), $(BLOCK), i) for i in xrange($(BLOCK), $(BLOCKCOUNT)+1)])')

BLKOUTFMT = $(BLOCKD)/$(BLOCKFMT)/$(BLOCKFMT).$(BLOCKFMT)
#DAFILEBASE := $(shell python -c 'print " ".join(["$(BLKOUTFMT)" % (i, j, i) for i in xrange(1,$(DABLOCKS)+1) for j in xrange(1, $(DABLOCKS)+1)])')
DAFILEBASE := $(shell python -c 'print " ".join(["$(BLKOUTFMT)" % ($(BLOCK), i, $(BLOCK)) for i in xrange(1, $(BLOCKCOUNT)+1)])')

DALPARTS = $(bp).N0.las $(bp).N1.las $(bp).N2.las $(bp).N3.las \
           $(bp).C0.las $(bp).C1.las $(bp).C2.las $(bp).C3.las
DALFILES :=  $(foreach bp, $(DAFILEBASE), $(DALPARTS))

# Use dazzler blocking for these steps
#BLOCKS := $(foreach b, $(shell seq 1 $(DABLOCKS)), $(BLOCKD)/$(SUBPFX).$(b))
BLOCKS := $(BLOCKD)/$(SUBPFX).$(BLOCK)
DZL1MRG := $(DAFILEBASE:%=%.las)
DZL2MRG := $(BLOCKS:%=%.las)
DAZ2M4 := $(BLOCKS:$(BLOCKD)/$(SUBPFX).%=correct/daz2m4.%.seeds)
MAPPEDM4 := $(DAZ2M4:correct/daz2m4.%.seeds=correct/seeds.%.m4)
MAPPEDM4FILT := $(MAPPEDM4:%=%.filt)
CORRECTED := $(MAPPEDM4:correct/seeds.%.m4=correct/corrected.%.fasta)

# Chunks are maintained throughout workflow. Need to add some elasticity to handle
# larger datasets
CHUNKS := $(shell seq 1 $(CHUNK_SIZE))
SPLITBESTN = $(shell echo $$(( 30/$(CHUNK_SIZE)+2 )))
BAXFOFNS := $(foreach c,$(CHUNKS),$(shell printf "input.chunk%03dof%03d.fofn" $(c) $(CHUNK_SIZE)))
REGFOFNS := $(BAXFOFNS:input.%=filter/regions.%)
CMPH5 := $(BAXFOFNS:input.%.fofn=polish/aligned_reads.%.cmp.h5)

CUTOFF = filter/longreads.cutoff

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
	$(QSUB) -N filt.$* "filterm4.py $< > $@"

# XXX: Fix this, really annoying, use fofn.  
# Because of dextract, minor edit to header and uppercase the bases, while joining into one file
filter/subreads.fasta : $(SUBFASTA)
	sed '/^[^>]/ y/acgt/ACGT/;s/ RQ=.*//' $^ > $@
##

## The key bridgepoint: adapt daligner output to m4 format ##
adapt : $(MAPPEDM4) ;

# Yeah, I wrote it in perl ... so what?
$(MAPPEDM4) : correct/seeds.%.m4 : $(SUBPFX).idmap $(BLOCKD)/$(SUBPFX).%.las
	$(QSUB) -N adpt.$* "daz2m4.pl $* $^ > $@"

##

## Merge dazzler alignments, arranged by target
dmerge : $(DZL2MRG) ;

# Modified LAMerge to remove limit of 252 files at a time.  Human typically generates more than that
$(DZL2MRG) : $(DZL1MRG)
	$(QSUB) -N d2m.$(@F) LAmerge $@ $(filter $(@:.las=/%),$^)

$(DZL1MRG) : %.las : %.N0.las %.N1.las %.N2.las %.N3.las %.C0.las %.C1.las %.C2.las %.C3.las | $(DALFILES)
	$(QSUB) -N d1s.$(@F) LAsort $^
	$(QSUB) -N d1m.$(@F) LAmerge $@ $(^:%.las=%.S.las)
	rm -f $(^:%.las=%.S.las)

#.INTERMEDIATE : $(DALFILES) $(DZL1MRG)
##

## Read overlap using DALIGNER
dalign : $(DALFILES) ;

# daligner will output both X.Y and Y.X, arrange block pairs by target for merging.
# XXX: Make sure the target resolves properly in the DAG after execution
$(DALFILES) : $(DALDONE)
	@if [ ! -e $@ ]; then ln -s $(PWD)/`echo $@ | sed -r 's;[0-9]+(/[a-z]+.)([0-9]+);\2\1\2;'` $@;else touch $@;fi

$(DALDONE) : correct/dbsplit.done | $(BLOCKS)
	cd $(@D) && $(QSUB) -N daln.$(@F) -pe smp 8 daligner -t7 `echo $(@F) | sed 's/_\([^ ]*\).done/ \1/'`
	touch $@

$(BLOCKS) : correct/$(SUBPFX).db | prepare
	mkdir -p $@
	ln -s $(PWD)/$< $@
	ln -s $(PWD)/correct/.$(SUBPFX).idx $@
	ln -s $(PWD)/correct/.$(SUBPFX).bps $@

##

## Find seed read cutoffs
$(CUTOFF) : $(SUBLENGTHS)
	sort -nrmk1,1 $^ | awk '{t+=$$1;if(t>=$(GENOME_SIZE)*30){print $$1;exit;}}' > $@

# NOTE: this also stores offsets that will eventually be used to map daz ids to pbids
$(SUBLENGTHS) : filter/%.$(SUBPFX).lengths : filter/%.$(SUBPFX).fasta
	fastalength $< | sort -nrk1,1 > $@

##

## Integration with dazzler flow, this has some explicit block enforcements
#  limited by the size of the daz internal read id.  All movies must be loaded
#  independently.  After creating the database, it must be split into blocks,
#  or chunks, if the # of reads is greater than uint16 (65,536) or number of
#  of bases is greater than some value (default 400Mbp).  Note This block
#  scheme will not be the same as the HGAP chunk scheme, so we'll have to 
#  manage dazzler blocks with HGAP chunks.
dazdb : $(SUBPFX).idmap correct/dbsplit.done

$(SUBPFX).idmap : $(CUTOFF) correct/$(SUBPFX).db 
	DB2idmap -c`cat $<` $(lastword $^)

correct/dbsplit.done : correct/$(SUBPFX).db
	DBsplit $< && touch $@

correct/$(SUBPFX).db : $(SUBFASTA)
	$(QSUB) -N f2db fasta2DB $@ $^

subreads : $(SUBFASTA) ;

$(SUBFASTA) : filter/%.$(SUBPFX).fasta : %.1.bax.h5 %.2.bax.h5 %.3.bax.h5 | prepare
	$(QSUB) -N dext.$* "dextract -s800 $^ > $@"

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

$(BAXFILES) : ;

##

clean :
	rm -rf $(TASKDIRS)
	rm -f input.chunk*
	rm -f rechunk*.done
	rm -f chunkinput.done
