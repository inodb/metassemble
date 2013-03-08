# Input parameters
ASSEMBLY_MAKEFILE=/bubo/home/h16/inod/glob/github/metassemble/scripts/Makefile
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
RUNNUCMERRULE=$(call get_sbatch_job_par,nucmer,-p node -t 02:00:00) \
	    "mkdir -p $(@D); \
		    bash $(SCRIPTDIR)/validate/nucmer/run-nucmer.sh $(REF) $< $(@D)/nucmer"
%minimus2/val/nucmer.coords: %minimus2/all-merged.fasta $(REF)
	sbatch $(RUNNUCMERRULE)
%newbler/val/nucmer.coords: %newbler/454AllContigs.fna $(REF)
	sbatch $(RUNNUCMERRULE)
%bambus2/val/nucmer.coords: %bambus2/bambus2.scaffold.linear.fasta $(REF)
	sbatch $(RUNNUCMERRULE)
%/val/nucmer.coords: %/contigs.fa $(REF)
	sbatch $(RUNNUCMERRULE)
%/val/nucmer.coords: %/meta-velvetg.contigs.fa $(REF)
	sbatch $(RUNNUCMERRULE)
%/valscaf/nucmer.coords: %/Scaffolds.fasta $(REF)
	sbatch $(RUNNUCMERRULE)
%/val/nucmer.coords: %/Contigs.fasta $(REF)
	sbatch $(RUNNUCMERRULE)

# Only validate existing, often something goes wrong in the assembly process so
# this allows you to only validate those assemblies that succeeded in case one
# is not able to get the other assemblies to complete.
.PHONY:
validateexisting: \
	$(subst bambus2.scaffold.linear.fasta,val/nucmer.coords,\
	$(subst Contigs.fasta,val/nucmer.coords,\
	$(subst Scaffolds.fasta,valscaf/nucmer.coords,\
	$(subst contigs.fa,val/nucmer.coords,\
	$(subst 454AllContigs.fna,val/nucmer.coords,\
	$(subst all-merged.fasta,val/nucmer.coords,\
	$(subst meta-velvetg.contigs.fa,val/nucmer.coords,$(wildcard $(ALLASMCONTIGS) $(ALLASMSCAFFOLDS)))))))))

.PRECIOUS: %minimus2/val/nucmer.coords %newbler/val/nucmer.coords %bambus2/val/nucmer.coords %/val/nucmer.coords
