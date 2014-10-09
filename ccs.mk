# Makefile for easy CCS shepherding
# Run using something like "make -j -f ccs.mk INPUT=path/input.fofn"
# dalexander, 2014

SHELL                   := /bin/bash
PWD                      = $(shell pwd)
QUEUE                    = primary

#
# Options
#
OUTPUT                  ?= MyCCS
INPUT                   ?= input.fofn
CCS_OPTIONS             ?= -v --minPredictedAccuracy=0
CHEMISTRY_OVERRIDE      ?=
OUTPUT_CSV              ?= --csv
REFERENCE               ?=
PARAMETERS_FILE         ?=

#
# Commands
#
QSUB                     = qsub -cwd -b y -sync y -V -q $(QUEUE) -e log/ -o log/
QSUB.1                   = $(QSUB) -pe smp 1
QSUB.4                   = $(QSUB) -pe smp 4
QSUB.8                   = $(QSUB) -pe smp 8

CCS                      = $(QSUB.8) ConsensusTools.sh CircularConsensus -n 8 $(CCS_OPTIONS) $(PARAMETERS_FILE) $(CHEMISTRY_OVERRIDE) $(OUTPUT_CSV)
MAP                      = $(QSUB.8) pbalign --nproc 8
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
MAP_OUTPUT              := $(CCS_OUTPUT:$(OUTPUT)/%.ccs.h5=$(OUTPUT)/%.cmp.h5)

$(OUTPUT)/%.ccs.h5 : %.bax.h5
	$(CCS) -o $(OUTPUT) $<

$(OUTPUT)/%.cmp.h5 : $(OUTPUT)/%.ccs.h5
	$(MAP_CCS) $< $(REFERENCE) $@

ccs: $(CCS_OUTPUT)

map: $(MAP_OUTPUT)

all: map

.PHONY: map ccs

