ifndef METASSEMBLE_DIR
$(error METASSEMBLE_DIR environment variable not set. Set with export METASSEMBLE_DIR=...)
endif
include $(METASSEMBLE_DIR)/scripts/parameters.mk
include $(METASSEMBLE_DIR)/lib/scheduler.mk

################################
# ----- general rules -------- #
################################
%.fastq.gz: %.fastq
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_GUNZIP_OPT),make -e $@)
################################
# ---- /general rules -------- #
################################

################################
# ------ quality trim -------- #
################################
$(FASTQ_TRIM_IL): $(FASTQ1) $(FASTQ2)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_QTRIM_OPT),make -e $@)
################################
# ------ /quality trim ------- #
################################

################################
# ---------- velveth --------- #
################################
$(VELVETH_OUT)/velveth_$(KMIN)/Sequences: $(FASTQ_TRIM_IL)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_VELVETH_OPT),make -e $@)
$(VELVETH_OUT)/velveth_%/Roadmaps: $(VELVETH_OUT)/velveth_$(KMIN)/Sequences
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_VELVETH_OPT),make -e $@)
################################
# --------- /velveth --------- #
################################

################################
# --------- velvetg ---------- #
################################
$(VELVET_OUT_NOSCAF)/noscaf_%/$(CONTIG_FILENAME): $(VELVETH_OUT)/velveth_%/Roadmaps
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_VELVETG_OPT),make -e $@)
$(VELVET_OUT_SCAF)/scaf_%/$(SCAF_FILENAME): $(VELVETH_OUT)/velveth_%/Roadmaps
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_VELVETG_OPT),make -e $@)
################################
# --------- /velvetg --------- #
################################

################################
# ------- meta-velvetg ------- #
################################
$(METAVELVET_OUT_NOSCAF)/noscaf_%/$(CONTIG_FILENAME): $(VELVETH_OUT)/velveth_%/Roadmaps
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_METAVELVETG_OPT),make -e $@)
$(METAVELVET_OUT_SCAF)/scaf_%/$(SCAF_FILENAME): $(VELVETH_OUT)/velveth_%/Roadmaps
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_METAVELVETG_OPT),make -e $@)
################################
# ------- /meta-velvetg ------ #
################################

################################
# --------- minimus2  -------- #
################################
$(MINIMUS2_OUT_VELVET_NOSCAF)/$(MERGE_FILENAME): $(VELVETG_OUT_NOSCAF)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_MINIMUS2_VELVET_OPT),make -e $@)
$(MINIMUS2_OUT_METAVELVET_NOSCAF)/$(MERGE_FILENAME): $(METAVELVETG_OUT_NOSCAF)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_MINIMUS2_METAVELVET_OPT),make -e $@)
$(MINIMUS2_OUT_RAY_NOSCAF)/$(MERGE_FILENAME): $(RAY_CONTIGS_OUT)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_MINIMUS2_RAY_OPT),make -e $@)
################################
# --------- /minimus2  ------- #
################################

################################
# --------- newbler -----------#
################################
$(NEWBLER_OUT_VELVET_NOSCAF)/$(MERGE_FILENAME): $(VELVETG_OUT_NOSCAF)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_NEWBLER_VELVET_OPT),make -e $@)
$(NEWBLER_OUT_METAVELVET_NOSCAF)/$(MERGE_FILENAME): $(METAVELVETG_OUT_NOSCAF)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_NEWBLER_METAVELVET_OPT),make -e $@)
$(NEWBLER_OUT_RAY_NOSCAF)/$(MERGE_FILENAME): $(RAY_CONTIGS_OUT)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_NEWBLER_RAY_OPT),make -e $@)
################################
# -------- /newbler -----------#
################################

################################
# --------- bambus2 -----------#
################################
%/bambus2/bambus2.scaffold.linear.fasta: %/$(MERGE_FILENAME)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_BAMBUS2_OPT),make -e $@)
%/bambus2/bambus2.scaffold.linear.fasta: %/$(CONTIG_FILENAME)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_BAMBUS2_OPT),make -e $@)
################################
# -------- /bambus2 -----------#
################################

################################
# ----------- ray -------------#
################################
$(RAY_OUT)/noscaf/noscaf_%/$(CONTIG_FILENAME): $(FASTQ_TRIM_IL)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_RAY_OPT),make -e $@)
$(RAY_OUT)/scaf/scaf_%/$(SCAF_FILENAME): $(RAY_OUT)/noscaf/noscaf_%/$(CONTIG_FILENAME)
	$(call schedule_with_deps_and_store_id,$(SCHEDULER_STD_OPT) $(SCHEDULER_RAY_MV_SCAF_OPT),make -e $@)
################################
# ---------- /ray -------------#
################################
