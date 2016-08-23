# A super basic implementation of HGAP.3 adapted for SEQUEL data. This 
# requires access to different major versions of smrtanalysis (2.3.0 and 3.x).
# You must source a 2.3.0 smrtanalysis build prior to running this script.  
# You must also configure a path to you 3.x.x smrtcmds bin directory.  See the 
# `bin3xx' setting below.
#
# It is slow, since it doesn't have any process-level parallelization defined, 
# but is useable for bacterial-sized genomes.
# Threads are defined at individual target levels, change to suit your needs.

### Parameters for adjustment

#-- REQUIRED
# path to a subreadset.xml that defines location and access to a BAM file
xml ?= undefined 

# You must have smrtanalysis 2.3.0 sourced and in your path for this to work.
# Modify path to respective bin of smrtanalysis 3.x.x
bin3xx := /pbi/dept/secondary/builds/mainline/current_smrttools_prebuilt_installdir/smrtcmds/bin

#-- DOUBLE CHECK
# defaults to an ecoli-type genome size
genome_size ?= 5000000
# cluster submission
qsub = qsub -V -cwd -q production -j y -b y -sync y -N $@ -pe smp ${threads}

#-- ADVANCED, defaults probably OK. Don't change unless you know what you're doing.
blasr_opts = -m 4 -nproc ${threads} -bestn 24 -noSplitSubreads -minReadLength 200 -maxScore -1000 -maxLCPLength 16
smrtetc = ${SEYMOUR_HOME}/analysis/etc

### End Parameters

SHELL := /bin/bash

fname := $(lastword ${MAKFILE_LIST})
subreads := subreads.fasta

all: polished.fasta

${xml}: 
	@echo "USAGE: make -f ${fname} xml=<path to subreadset.xml>" && exit 1

# adapt from sequel format
${subreads}.gz: ${xml}
	${bin3xx}/bam2fasta -o $(basename ${subreads}) $<

${subreads}: subreads.fasta.gz
	gunzip -c $< > $@

# get seed read cutoff
length.cutoff: ${subreads}
	fastalength $< | sort -nrmk1,1 | awk '{t+=$$1;if(t>=${genome_size}*30){print $$1;exit;}}' > $@
	
# partition seed (long) reads
longreads.fasta: length.cutoff
	fastalength $< | awk -v len=$$(cat $<) '($$1<len){print $$2}' | fastaremove ${subreads} stdin > $@

# overlap seeds against all reads
overlaps.m4: threads = 24
overlaps.m4: longreads.fasta
	${qsub} blasr ${subreads} $< ${blasr_opts} -out $@

# filter overlaps for chimeras, artifacts and such
overlaps.filt.m4: overlaps.m4
	filterm4.py $< > $@

# correct seeds via consensus
corrected.fa: threads = 15
corrected.fa: overlaps.filt.m4
	${qsub} -v mym4=$<,allm4=$<,subreads=${subreads},cov=6,tmp=/scratch pbdagcon_wf.sh

corrected.fq: corrected.fa

# assemble corrected seeds into draft
utg.spec: corrected.fa
	runCASpecWriter.py  -vv --bitTable=${smrtetc}/celeraAssembler/bitTable \
	--interactiveTmpl=${smrtetc}/cluster/SGE/interactive.tmpl \
	--smrtpipeRc=${smrtetc}/smrtpipe.rc --genomeSize=5000000 --defaultFrgMinLen=500 \
	--xCoverage=20 --ovlErrorRate=0.06 --ovlMinLen=40 --merSize=14 --corrReadsFasta=$< \
	--specOut=$@ --sgeName=utg --gridParams="useGrid:0, scriptOnGrid:0, frgCorrOnGrid:0, ovlCorrOnGrid:0" \
	--maxSlotPerc=1 ${smrtetc}/celeraAssembler/unitig.spec

utg.frg: corrected.fq
	fastqToCA -technology sanger -type sanger -libraryname corr -reads $< > $@

utg.finished: threads = 16
utg.finished: utg.spec utg.frg
	${qsub} -pe smp 16 runCA -d assemble -p asm -s $^
	touch $@

draft.fasta: threads = 4
draft.fasta: utg.finished
	tigStore -g assemble/asm.gkpStore -t assemble/asm.tigStore 1 -d properties -U \
	| awk 'BEGIN{t=0}$$1=="numFrags"{if($$2>1){print t, $$2}t++}' | sort -nrk2,2 > assemble/unitig.lst
	${qsub} 'tmp=$$(mktemp -d -p /scratch); tmp=$$tmp cap=${PWD}/assemble/asm utg=${PWD}/assemble/unitig.lst cns=${PWD}/$@ nproc=${threads} pbutgcns_wf.sh'

# prepare for polishing
draft.fasta.fai: draft.fasta
	samtools faidx $<

draft.mapped.bam: threads = 24
draft.mapped.bam: ${xml} draft.fasta
	${qsub} ${bin3xx}/blasr $^ --nproc ${threads} --bestn 1 --maxScore -1000 --hitPolicy randombest --bam --out $@

draft.sorted.bam: draft.mapped.bam
	samtools sort $< $(basename $@)

draft.sorted.bam.pbi: draft.sorted.bam
	${bin3xx}/pbindex $<

# polish
polished.fasta: threads = 16
polished.fasta: draft.fasta draft.sorted.bam draft.fasta.fai draft.sorted.bam.pbi
	${qsub} ${bin3xx}/arrow -j ${threads} -o $@ --referenceFilename $(wordlist 1,2,$^)

clean:
	rm -f subreads.fasta	
	rm -f subreads.fasta.gz
	rm -f length.cutoff
	rm -f longreads.fasta
	rm -f corrected.f*
	rm -f draft.*
	rm -f overlaps.*
	rm -f utg.*
	rm -rf assemble/
	rm -f polished.fasta