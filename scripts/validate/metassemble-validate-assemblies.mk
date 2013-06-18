# Input parameters
ASSEMBLY_MAKEFILE=$(METASSEMBLE_DIR)/scripts/metassemble.mk
REF=/bubo/home/h16/inod/metagenomics/results/chris-mock/Project_ID793_dAmore/Sample_50ng_even/mapping/velvet35test/references.fa

include $(ASSEMBLY_MAKEFILE)

# Calculate reference stats
$(OUT)/reference-stats/ref.stats:
	sbatch $(call get_sbatch_job_par,mapref,-p node -t 2-00:00:00) \
		"mkdir -p $(dir $@); \
		bash $(SCRIPTDIR)/validate/reference/stats/length-gc-cov.sh \
			$(PRC_READS_OUT)/$(FASTQBASE).1.qtrim $(PRC_READS_OUT)/$(FASTQBASE).2.qtrim \
			$(FASTQBASE).qtrim $(REF) ref $(dir $@) > $@"

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

%/val/asm-stats.tsv: %/bambus2.scaffold.linear.fasta $(REF)
	$(MASMVALIRULE)
%/val/asm-stats.tsv: %/$(CONTIG_FILENAME) $(REF)
	$(MASMVALIRULE)
%/val/asm-stats.tsv: %/$(SCAF_FILENAME) $(REF)
	$(MASMVALIRULE)
%/val/asm-stats.tsv: %/$(MERGE_FILENAME) $(REF)
	$(MASMVALIRULE)

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
