#
# Makefile for "titrating" different coverage levels through quiver
# caller to check convergence of consensus accuracy.
#
# Run using something like ""
#
# [Awaiting easier readscore filtering in pbalign before providing support
#  for starting from bax files]
#
# Present legal values of GENOME:
#   - ecoliK12_pbi_March2013
#   - S_aureus_USA300_TCH1516
#   - R_palustris_CGA009
#
# dalexander, 2014

SHELL                   := /bin/bash
PWD                      = $(shell pwd)
QUEUE                    = secondary


# ------------------------------------------------------------------------------
# Options
#
CMPH5                   ?=
#INPUT_FOFN              ?=
CHEMISTRY_OVERRIDE      ?=
GENOME                  ?=
CONDITION               ?= $(GENOME)
COVERAGE_LEVELS         ?= 5 10 15 20 30 40 60 80
PARAMETERS_FILE         ?=


ifndef GENOME
$(error GENOME is not set)
endif

ifndef CMPH5
$(error CMPH5 is not set)
endif

# Root paths for:
# 1) the pacbio "reference repository" where reference FASTAs live
# 2) a repository of GFF files indicating portions of the reference
#    FASTAs that are known to be inaccurate and so should be excluded
#    from analysis
# These two will be keyed by the $(GENOME) variable, so it is important
# to make sure that they match up!
REFERENCE_REPOSITORY_ROOT := /mnt/secondary/iSmrtanalysis/current/common/references
REFERENCE_MASK_ROOT       := /mnt/secondary/Share/Quiver/GenomeMasks


# TODO: covariates??


# ------------------------------------------------------------------------------
# Computed variables
#
REFERENCE               := $(REFERENCE_REPOSITORY_ROOT)/$(GENOME)/sequence/$(GENOME).fasta
REFERENCE_MASK          := $(REFERENCE_MASK_ROOT)/$(GENOME)-mask.gff
OUTPUT                  ?= Output/$(CONDITION)

RAW_VARIANTS            := $(foreach COVERAGE,$(COVERAGE_LEVELS),$(OUTPUT)/RawVariants/$(CONDITION).$(COVERAGE).gff)
MASKED_VARIANTS         := $(foreach COVERAGE,$(COVERAGE_LEVELS),$(OUTPUT)/MaskedVariants/$(CONDITION).$(COVERAGE).gff)


# ------------------------------------------------------------------------------
# Commands
#
QSUB                     = qsub -cwd -b y -sync y -V -q $(QUEUE) -e log/ -o log/
QSUB.1                   = $(QSUB) -pe smp 1
QSUB.4                   = $(QSUB) -pe smp 4
QSUB.8                   = $(QSUB) -pe smp 8
QSUB.16                  = $(QSUB) -pe smp 16

# Turn of all variant filtering so we see all discrepancies from the
# reference.  TODO: output nocalls as "variant" type under flag?
QUIVER                   = $(QSUB.8) quiver -j8 -v -x0 -q0  $(PARAMETERS_FILE) $(CHEMISTRY_OVERRIDE)


# ------------------------------------------------------------------------------
# Setup
MKDIR_OK                := $(shell mkdir -p $(OUTPUT) && \
	                     mkdir -p $(OUTPUT)/RawVariants && \
	                     mkdir -p $(OUTPUT)/MaskedVariants && \
			     mkdir -p log)

# ------------------------------------------------------------------------------
# Targets

coverage = $(shell echo $(1) | awk -F. '{print $$(NF-1)}')

titration: $(OUTPUT)/titration.csv

$(OUTPUT)/titration.csv: $(MASKED_VARIANTS)
	scripts/makeTitrationCsv.sh

$(OUTPUT)/MaskedVariants/%.gff : $(OUTPUT)/RawVariants/%.gff $(REFERENCE_MASK)
	gffsubtract.pl $< $(REFERENCE_MASK) > $@

$(OUTPUT)/RawVariants/%.gff : $(CMPH5) $(REFERENCE)
	$(QUIVER) -X $(call coverage,$@) $(CMPH5) -r $(REFERENCE) -o $@

.PHONY: titration map

.PRECIOUS: $(RAW_VARIANTS)
