FASTA ?= polished_assembly.fasta
OVLSIZE = 20000
ODIR = circularize

all :  circ

prepare :
	@mkdir -p $(ODIR)

circ : $(ODIR)/ovl.m4
	@awk '($$6==0 && $$9==0 && $$12==$(OVLSIZE)){print substr($$1,0,index($$1,":"))" looks circular"}' $<

$(ODIR)/ovl.m4 : $(ODIR)/beg.fa $(ODIR)/end.fa
	@blasr $^ -bestn 1 -m 4 -minPctIdentity 98.0 > $@ 2> /dev/null

$(ODIR)/beg.fa : $(FASTA) | prepare
	@fastasubseq $< 0 $(OVLSIZE) > $@

$(ODIR)/end.fa : $(FASTA) $(ODIR)/end.offset | prepare
	@fastasubseq $< `cat $(word 2,$^)` $(OVLSIZE) > $@

$(ODIR)/end.offset : $(FASTA) | prepare
	@fastalength $< | awk '{print $$1-$(OVLSIZE)}' > $@

clean :
	@rm -rf $(ODIR)
