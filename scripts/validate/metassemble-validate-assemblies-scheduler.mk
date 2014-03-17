ifndef METASSEMBLE_DIR
$(error METASSEMBLE_DIR environment variable not set. Set with export METASSEMBLE_DIR=...)
endif
include $(METASSEMBLE_DIR)/scripts/metassemble.mk
include $(METASSEMBLE_DIR)/lib/scheduler.mk

ifndef REF
$(error REF variable not set. Set like REF=/path/to/references.fa)
endif
ifndef PHYL_REF
$(error PHYL_REF variable not set. Set like PHYL_REF=/path/to/references.fa)
endif
ifndef MAKEFILE_VALIDATE
$(error MAKEFILE_VALIDATE variable not set. Set like MAKEFILE_VALIDATE=/path/to/references.fa)
endif

include $(ASSEMBLY_MAKEFILE)

################################
# --------- ref stats ---------#
################################
$(OUT)/reference-stats/ref.stats: $(FASTQ1) $(FASTQ2) $(REF)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_REFMAP_OPT),make -ef $(MAKEFILE_VALIDATE) $@)
################################
# ---------/ref stats ---------#
################################

################################
# --------- nucmer ------------#
################################
%/val/nucmer.coords: %/bambus2.scaffold.linear.fasta $(REF)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_NUCMER_OPT),make -ef $(MAKEFILE_VALIDATE) $@)
%/val/nucmer.coords: %/$(CONTIG_FILENAME) $(REF)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_NUCMER_OPT),make -ef $(MAKEFILE_VALIDATE) $@)
%/val/nucmer.coords: %/$(SCAF_FILENAME) $(REF)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_NUCMER_OPT),make -ef $(MAKEFILE_VALIDATE) $@)
%/val/nucmer.coords: %/$(MERGE_FILENAME) $(REF)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_NUCMER_OPT),make -ef $(MAKEFILE_VALIDATE) $@)

%/val/asm-stats.tsv: %/bambus2.scaffold.linear.fasta %/val/nucmer.coords $(REF) $(PHYL_REF) $(CON_TO_REF) $(OUT)/reference-stats/ref.stats
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_MASMVALI_OPT),make -ef $(MAKEFILE_VALIDATE) $@)
%/val/asm-stats.tsv: %/$(CONTIG_FILENAME) %/val/nucmer.coords $(REF) $(PHYL_REF) $(CON_TO_REF) $(OUT)/reference-stats/ref.stats
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_MASMVALI_OPT),make -ef $(MAKEFILE_VALIDATE) $@)
%/val/asm-stats.tsv: %/$(SCAF_FILENAME) %/val/nucmer.coords $(REF) $(PHYL_REF) $(CON_TO_REF) $(OUT)/reference-stats/ref.stats
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_MASMVALI_OPT),make -ef $(MAKEFILE_VALIDATE) $@)
%/val/asm-stats.tsv: %/$(MERGE_FILENAME) %/val/nucmer.coords $(REF) $(PHYL_REF) $(CON_TO_REF) $(OUT)/reference-stats/ref.stats
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_MASMVALI_OPT),make -ef $(MAKEFILE_VALIDATE) $@)
################################
# ---------/nucmer ------------#
################################

################################
# --------- only map ----------#
################################
%/val/map/bowtie2/asm_$(FASTQBASE)-smds.bam: %/$(CONTIG_FILENAME) $(FASTQ_TRIM_1) $(FASTQ_TRIM_2)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_MAP_BOWTIE_OPT),make -ef $(MAKEFILE_VALIDATE) $@)
%/val/map/bowtie2/asm_$(FASTQBASE)-smds.bam: %/$(MERGE_FILENAME) $(FASTQ_TRIM_1) $(FASTQ_TRIM_2)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_MAP_BOWTIE_OPT),make -ef $(MAKEFILE_VALIDATE) $@)
%/val/map/bowtie2/asm_$(FASTQBASE)-smds.bam: %/$(SCAF_FILENAME) $(FASTQ_TRIM_1) $(FASTQ_TRIM_2)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_MAP_BOWTIE_OPT),make -ef $(MAKEFILE_VALIDATE) $@)
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
