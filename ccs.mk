# Makefile for easy CCS shepherding
# Run using something like "make -j -f ccs.mk INPUT=path/input.fofn"
# dalexander, 2014

SHELL                   := /bin/bash
PWD                      = $(shell pwd)
QUEUE                    =

#
# Options
#
OUTPUT                  ?= Output
INPUT                   ?= input.fofn
CCS_OPTIONS             ?= -v --minPredictedAccuracy=0
CHEMISTRY_OVERRIDE      ?=
OUTPUT_CSV              ?= --csv
REFERENCE               ?=
PARAMETERS_FILE         ?=

#
# Commands
#
QSUB                     = qsub -cwd -b y -sync y -V -e log/ -o log/
QSUB.1                   = $(QSUB) -pe smp 1
QSUB.4                   = $(QSUB) -pe smp 4
QSUB.8                   = $(QSUB) -pe smp 8

CCS                      = $(QSUB.8) ConsensusTools.sh CircularConsensus -n 8 $(CCS_OPTIONS) $(PARAMETERS_FILE) $(CHEMISTRY_OVERRIDE) $(OUTPUT_CSV)


MAPPING_OPTIONS          =
MAP                      = $(QSUB.8) pbalign --nproc 8 $(MAPPING_OPTIONS)
MAP_CCS                  = $(MAP) --useccs=useccsdenovo

#
# Setup before rule execution.  Note we make an assumption that the
# bax file "basename" is unique among all the directories in question,
# because of the VPATH approach.  Fix this in the future.
#
MKDIR_OK                := $(shell mkdir -p $(OUTPUT) && mkdir -p log)
VPATH                   := $(shell cat $(INPUT) | xargs -n1 dirname | sort | uniq)
BAX_FILES               := $(shell cat $(INPUT) | xargs -n1 basename | sort | uniq)

#
# Rules
#
CCS_OUTPUT              := $(BAX_FILES:%.bax.h5=$(OUTPUT)/%.ccs.h5)
MAP_OUTPUT              := $(BAX_FILES:%.bax.h5=$(OUTPUT)/%.subreads.cmp.h5)
MAP_CCS_OUTPUT          := $(CCS_OUTPUT:$(OUTPUT)/%.ccs.h5=$(OUTPUT)/%.ccs.cmp.h5)

$(OUTPUT)/%.ccs.h5 : %.bax.h5
	$(CCS) -o $(OUTPUT) $<

$(OUTPUT)/%.subreads.cmp.h5 : %.bax.h5
	$(MAP) $< $(REFERENCE) $@

$(OUTPUT)/%.ccs.cmp.h5 : $(OUTPUT)/%.ccs.h5
	$(MAP_CCS) $< $(REFERENCE) $@

ccs: $(CCS_OUTPUT)

map: $(MAP_OUTPUT)

map-ccs: $(MAP_CCS_OUTPUT)

all: map map-ccs

.PHONY: map ccs map-ccs all
