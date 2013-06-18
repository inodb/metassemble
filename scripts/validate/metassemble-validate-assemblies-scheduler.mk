REF?=/bubo/home/h16/inod/glob/metagenomics/result/chris-mock/Project_ID793_dAmore/Sample_50ng_even/mapping/velvet35test/references.fa
PHYL_REF?=/bubo/home/h16/inod/glob/github/masm-vali/masmvali/test/data/chris-mock/reference/phylogeny-references.tsv
MAKEFILE_VALIDATE?=Makefile-validate-assemblies

ifndef METASSEMBLE_DIR
$(error METASSEMBLE_DIR environment variable not set. Set with export METASSEMBLE_DIR=...)
endif
include $(METASSEMBLE_DIR)/scripts/metassemble.mk
include $(METASSEMBLE_DIR)/lib/scheduler.mk

%/val/nucmer.coords: %/$(CONTIG_FILENAME) $(REF)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_NUCMER_OPT),make -ef $(MAKEFILE_VALIDATE) $@)
%/val/nucmer.coords: %/$(SCAF_FILENAME) $(REF)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_NUCMER_OPT),make -ef $(MAKEFILE_VALIDATE) $@)
%/val/nucmer.coords: %/$(MERGE_FILENAME) $(REF)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_NUCMER_OPT),make -ef $(MAKEFILE_VALIDATE) $@)

%/val/asm-stats.tsv: %/val/nucmer.coords $(OUT)/reference-stats/ref.stats $(PHYL_REF) %/$(CONTIG_FILENAME)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_MASMVALI_OPT),make -ef $(MAKEFILE_VALIDATE) $@)
%/val/asm-stats.tsv: %/val/nucmer.coords $(OUT)/reference-stats/ref.stats $(PHYL_REF) %/$(SCAF_FILENAME)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_MASMVALI_OPT),make -ef $(MAKEFILE_VALIDATE) $@)
%/val/asm-stats.tsv: %/val/nucmer.coords $(OUT)/reference-stats/ref.stats $(PHYL_REF) %/$(MERGE_FILENAME)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_MASMVALI_OPT),make -ef $(MAKEFILE_VALIDATE) $@)

# Only validate existing, often something goes wrong in the assembly process so
# this allows you to only validate those assemblies that succeeded in case one
# is not able to get the other assemblies to complete.
.PHONY:
validateexisting: \
	$(subst $(MERGE_FILENAME),val/asm-stats.tsv,\
	$(subst $(SCAF_FILENAME),val/asm-stats.tsv,\
	$(subst $(CONTIG_FILENAME),val/asm-stats.tsv,$(wildcard $(ALLASMCONTIGS) $(ALLASMSCAFFOLDS)))))

.PRECIOUS: %/val/nucmer.coords
