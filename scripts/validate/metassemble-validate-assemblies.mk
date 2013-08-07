# Input parameters
ASSEMBLY_MAKEFILE=$(METASSEMBLE_DIR)/scripts/metassemble.mk
REF=/bubo/home/h16/inod/metagenomics/results/chris-mock/Project_ID793_dAmore/Sample_50ng_even/mapping/velvet35test/references.fa

include $(ASSEMBLY_MAKEFILE)


################################
# --------- ref stats ---------#
################################
# Calculate reference stats
$(OUT)/reference-stats/ref.stats: $(FASTQ1) $(FASTQ2) $(REF)
	mkdir -p $(dir $@)
	bash $(SCRIPTDIR)/validate/reference/stats/length-gc-cov.sh \
		$(FASTQ1) $(FASTQ2) \
		$(FASTQBASE) $(REF) ref $(dir $@) > $@
################################
# ---------/ref stats ---------#
################################


################################
# --------- nucmer ------------#
################################
# Validate assemblies, run nucmer
RUNNUCMERRULE=mkdir -p $(@D); \
bash $(SCRIPTDIR)/validate/nucmer/run-nucmer.sh $(REF) $< $(@D)/nucmer

%/val/nucmer.coords: %/bambus2.scaffold.linear.fasta $(REF)
	$(RUNNUCMERRULE)
%/val/nucmer.coords: %/$(CONTIG_FILENAME) $(REF)
	$(RUNNUCMERRULE)
%/val/nucmer.coords: %/$(SCAF_FILENAME) $(REF)
	$(RUNNUCMERRULE)
%/val/nucmer.coords: %/$(MERGE_FILENAME) $(REF)
	$(RUNNUCMERRULE)

# Validate assemblies, run nucmer
MASMVALIRULE=python /glob/inod/github/masm-vali/masmvali/validation.py $(@D)/nucmer.coords $(OUT)/reference-stats/ref.stats /glob/inod/github/masm-vali/masmvali/test/data/chris-mock/reference/phylogeny-references.tsv $< $(@D)

%/val/asm-stats.tsv: %/bambus2.scaffold.linear.fasta %/val/nucmer.coords $(REF)
	$(MASMVALIRULE)
%/val/asm-stats.tsv: %/$(CONTIG_FILENAME) %/val/nucmer.coords $(REF)
	$(MASMVALIRULE)
%/val/asm-stats.tsv: %/$(SCAF_FILENAME) %/val/nucmer.coords $(REF)
	$(MASMVALIRULE)
%/val/asm-stats.tsv: %/$(MERGE_FILENAME) %/val/nucmer.coords $(REF)
	$(MASMVALIRULE)
################################
# ---------/nucmer ------------#
################################

################################
# --------- only map ----------#
################################
define MAP_BOWTIE_RULE
mkdir -p $(@D)
bash -x $(SCRIPTDIR)/map/map-bowtie-markduplicates.sh $(MAP_PARS) -c $(wordlist 2, 2, $^) $(lastword $^) \
	pair $< asm $(@D)
endef	

%/val/map/bowtie2/asm_$(FASTQBASE)-smds.bam: %/$(CONTIG_FILENAME) $(FASTQ_TRIM_1) $(FASTQ_TRIM_2)
	$(MAP_BOWTIE_RULE)
%/val/map/bowtie2/asm_$(FASTQBASE)-smds.bam: %/$(MERGE_FILENAME) $(FASTQ_TRIM_1) $(FASTQ_TRIM_2)
	$(MAP_BOWTIE_RULE)
%/val/map/bowtie2/asm_$(FASTQBASE)-smds.bam: %/$(SCAF_FILENAME) $(FASTQ_TRIM_1) $(FASTQ_TRIM_2)
	$(MAP_BOWTIE_RULE)
################################
# ---------/only map ----------#
################################

# Only validate existing, often something goes wrong in the assembly process so
# this allows you to only validate those assemblies that succeeded in case one
# is not able to get the other assemblies to complete.
.PHONY:
validateexisting: \
	$(subst bambus2.scaffold.linear.fasta,val/asm-stats.tsv,\
	$(subst $(MERGE_FILENAME),val/asm-stats.tsv,\
	$(subst $(SCAF_FILENAME),val/asm-stats.tsv,\
	$(subst $(CONTIG_FILENAME),val/asm-stats.tsv,$(wildcard $(ALLASMCONTIGS) $(ALLASMSCAFFOLDS))))))

.PRECIOUS: %minimus2/val/nucmer.coords %newbler/val/nucmer.coords %bambus2/val/nucmer.coords %/val/nucmer.coords
