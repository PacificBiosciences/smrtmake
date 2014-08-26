FASTA ?= polished_assembly.fasta
OVLSIZE ?= 20000
ODIR = circularize
# sequence ids suitable for file names (i.e., parse out |quiver)
FIDS := $(shell fastalength $(FASTA) | awk '($$1>2*$(OVLSIZE)){print $$2}' | sed -r 's/([^|]+).*/\1/')

EOFF := $(FIDS:%=$(ODIR)/%.end.offs)
BEGS := $(FIDS:%=$(ODIR)/%.beg.fa)
ENDS := $(FIDS:%=$(ODIR)/%.end.fa)
FAS  := $(FIDS:%=$(ODIR)/%.fa)
M4S  := $(FIDS:%=$(ODIR)/%.m4)
CIRC := $(FIDS:%=%.circ)

all : $(CIRC)

prepare :
	@mkdir -p $(ODIR)

$(CIRC) : %.circ : $(ODIR)/%.m4 
	@awk '($$6==0 && $$9==0 && $$12==$(OVLSIZE)){print "$*: looks circular"}' $<

$(M4S) : %.m4 : %.beg.fa %.end.fa
	@blasr $^ -bestn 1 -m 4 -minPctIdentity 98.0 > $@ 2> /dev/null

$(BEGS) : %.beg.fa : %.fa
	@fastasubseq $< 0 $(OVLSIZE) > $@


$(ENDS) : %.end.fa : %.fa %.end.offs 
	@fastasubseq $< `cat $(word 2,$^)` $(OVLSIZE) > $@

$(EOFF) : %.end.offs : %.fa 
	@fastalength $< | awk '{print $$1-$(OVLSIZE)}' > $@

$(FAS) : $(ODIR)/%.fa : $(FASTA) % | prepare
	@awk '{if($$1~">"){p=$$1~/$*\y/?1:0}if(p){print}}' $< > $@

$(FIDS) : ;

clean :
	@rm -rf $(ODIR)

#.INTERMEDIATE : $(EOFF) $(FAS) $(BEGS) $(ENDS) $(M4S)
