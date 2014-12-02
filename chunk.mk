# Provides a simple recipe for creating chunks out of bax file input.fofn
# Include this in your makefile and use the *_CHUNKS variables to create more 
# recipes.

BAXFOFN ?= input.fofn
CHUNK ?= 3

CHUNK_FMT = chunk%03dof%03d
INIT_CHUNKS := $(foreach c,$(shell seq 1 $(CHUNK)),$(shell printf ".init.$(CHUNK_FMT)" $(c) $(CHUNK)))
BAXFOFN_CHUNKS := $(INIT_CHUNKS:.init.%=input.%.fofn)
N = $(shell echo $* | sed -r 's/chunk0?0?([0-9]+)of[0-9]+/\1/')

all : $(BAXFOFN_CHUNKS)

$(BAXFOFN_CHUNKS) : input.%.fofn : $(BAXFOFN) .init.%
	awk 'BEGIN{c=$(N)%$(CHUNK)}(NR%$(CHUNK)==c){print $$0}' $< > $@

$(INIT_CHUNKS) : $(BAXFOFN) 
	@touch $@

$(BAXFOFN) : ;

