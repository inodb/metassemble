# Input parameters
ifndef METASSEMBLE_DIR
$(error METASSEMBLE_DIR environment variable not set. Set with export METASSEMBLE_DIR=...)
endif
ASSEMBLY_MAKEFILE=$(METASSEMBLE_DIR)/scripts/metassemble.mk
ifndef REF
$(error REF variable not set. Set like REF=/path/to/references.fa)
endif
ifndef PHYL_REF
$(error PHYL_REF variable not set. Set like PHYL_REF=/path/to/references.fa)
endif
ifndef CON_TO_REF
$(error CON_TO_REF variable not set. Set like CON_TO_REF=/path/to/contigs_to_refs.tsv)
endif

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
MASMVALIRULE=python /glob/inod/github/masm-vali/masmvali/validation.py --contigs_to_refs_table $(CON_TO_REF) $(@D)/nucmer.coords $(OUT)/reference-stats/ref.stats $(PHYL_REF) $< $(@D)

%/val/asm-stats.tsv: %/bambus2.scaffold.linear.fasta %/val/nucmer.coords $(REF) $(PHYL_REF) $(CON_TO_REF) $(OUT)/reference-stats/ref.stats
	$(MASMVALIRULE)
%/val/asm-stats.tsv: %/$(CONTIG_FILENAME) %/val/nucmer.coords $(REF) $(PHYL_REF) $(CON_TO_REF) $(OUT)/reference-stats/ref.stats
	$(MASMVALIRULE)
%/val/asm-stats.tsv: %/$(SCAF_FILENAME) %/val/nucmer.coords $(REF) $(PHYL_REF) $(CON_TO_REF) $(OUT)/reference-stats/ref.stats
	$(MASMVALIRULE)
%/val/asm-stats.tsv: %/$(MERGE_FILENAME) %/val/nucmer.coords $(REF) $(PHYL_REF) $(CON_TO_REF) $(OUT)/reference-stats/ref.stats
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
.PHONY:
validateall: \
	$(subst bambus2.scaffold.linear.fasta,val/asm-stats.tsv,\
	$(subst $(MERGE_FILENAME),val/asm-stats.tsv,\
	$(subst $(SCAF_FILENAME),val/asm-stats.tsv,\
	$(subst $(CONTIG_FILENAME),val/asm-stats.tsv,$(ALLASMCONTIGS) $(ALLASMSCAFFOLDS)))))

.PRECIOUS: %/val/nucmer.coords
