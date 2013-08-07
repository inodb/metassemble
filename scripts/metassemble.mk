# MetAssemble Pipeline. This Makefile is NOT meant to build the scripts but to
# run the metassemble scripts on an Illumina paired end library. An example of
# usage can be found in the examples/Makefile.

ifndef METASSEMBLE_DIR
$(error METASSEMBLE_DIR environment variable not set. Set with export METASSEMBLE_DIR=...)
endif
include $(METASSEMBLE_DIR)/scripts/parameters.mk
SCRIPTDIR=$(METASSEMBLE_DIR)/scripts

################################
# ----- general rules -------- #
################################
%.fastq.gz: %.fastq
	gzip $<
################################
# ---- /general rules -------- #
################################

################################
# ------ quality trim -------- #
################################
ifeq ($(DO_QTRIM),yes)
$(FASTQ_TRIM_1): $(FASTQ1) $(FASTQ2)
	mkdir -p $(PRC_READS_OUT)
	sickle pe \
		-f $(FASTQ1) \
		-r $(FASTQ2) \
		-t sanger \
		-o $(FASTQ_TRIM_1) \
		-p $(FASTQ_TRIM_2) \
		-s $(FASTQ_TRIM_IL).unpaired
$(FASTQ_TRIM_2): $(FASTQ_TRIM_1)
	if [ -f $@ ]; then \
		touch $@; \
	else \
		exit 1; \
	fi
endif
$(FASTQ_TRIM_IL): $(FASTQ_TRIM_1) $(FASTQ_TRIM_2)
	mkdir -p $(PRC_READS_OUT)
	shuffleSequences_fastq.pl $(FASTQ_TRIM_1) \
		$(FASTQ_TRIM_2) \
		$@
################################
# ------ /quality trim ------- #
################################

################################
# ---------- velveth --------- #
################################
$(VELVETH_OUT)/velveth_$(KMIN)/Sequences: $(FASTQ_TRIM_IL)
	mkdir -p $(VELVETH_OUT)
	velveth $(VELVETH_OUT)/velveth_$(KMIN) $(KMIN) -noHash -fastq \
		-shortPaired $<
$(VELVETH_OUT)/velveth_%/Sequences: $(VELVETH_OUT)/velveth_$(KMIN)/Sequences
	mkdir -p $(@D)
	ln -fs $(abspath $<) $@
$(VELVETH_OUT)/velveth_%/Roadmaps: $(VELVETH_OUT)/velveth_%/Sequences
	velveth $(@D) $* -reuse_Sequences
################################
# --------- /velveth --------- #
################################

################################
# --------- velvetg ---------- #
################################
define velvetg_rule
mkdir -p $(@D)
ln -fs $(abspath $(lastword $^)) $(@D)/Sequences
ln -fs $(abspath $<) $(@D)/Roadmaps
velvetg $(@D) $1
mv $(@D)/contigs.fa $@
endef
$(VELVET_OUT_NOSCAF)/noscaf_%/$(CONTIG_FILENAME): $(VELVETH_OUT)/velveth_%/Roadmaps $(VELVETH_OUT)/velveth_%/Sequences
	$(call velvetg_rule,-scaffolding no)
$(VELVET_OUT_SCAF)/scaf_%/$(SCAF_FILENAME): $(VELVETH_OUT)/velveth_%/Roadmaps $(VELVETH_OUT)/velveth_%/Sequences
	$(call velvetg_rule,-scaffolding yes -exp_cov auto)
################################
# --------- /velvetg --------- #
################################

################################
# ------- meta-velvetg ------- #
################################
# Copy output from velveth and run velvetg, followed by meta-velvetg -scaffolding yes or no
define metavelvetg_rule
mkdir -p $(dir $@)
ln -fs $(abspath $(lastword $^)) $(@D)/Sequences
ln -fs $(abspath $<) $(@D)/Roadmaps
velvetg $(dir $@) -scaffolding no -exp_cov auto -read_trkg yes \
	&& meta-velvetg $(dir $@) $1
mv $(@D)/meta-velvetg.contigs.fa $@
endef
$(METAVELVET_OUT_NOSCAF)/noscaf_%/$(CONTIG_FILENAME): $(VELVETH_OUT)/velveth_%/Roadmaps $(VELVETH_OUT)/velveth_%/Sequences
	$(call metavelvetg_rule,-scaffolding no)
$(METAVELVET_OUT_SCAF)/scaf_%/$(SCAF_FILENAME): $(VELVETH_OUT)/velveth_%/Roadmaps $(VELVETH_OUT)/velveth_%/Sequences
	$(call metavelvetg_rule,-scaffolding yes)
################################
# ------- /meta-velvetg ------ #
################################

################################
# --------- minimus2  -------- #
################################
# Minimus2 rule merges all given prerequisite files
define MINIMUS2_RULE
mkdir -p $(@D)
bash -x $(SCRIPTDIR)/assembly/merge-asm-minimus2.sh $(@D) $^
mv $(@D)/all-merged.fasta $@
endef
$(MINIMUS2_OUT_VELVET_NOSCAF)/$(MERGE_FILENAME): $(VELVETG_OUT_NOSCAF)
	$(MINIMUS2_RULE)
$(MINIMUS2_OUT_METAVELVET_NOSCAF)/$(MERGE_FILENAME): $(METAVELVETG_OUT_NOSCAF)
	$(MINIMUS2_RULE)
$(MINIMUS2_OUT_RAY_NOSCAF)/$(MERGE_FILENAME): $(RAY_CONTIGS_OUT)
	$(MINIMUS2_RULE)
################################
# --------- /minimus2  ------- #
################################

################################
# --------- newbler -----------#
################################
# Newbler rule merges all given prerequisite files
define NEWBLER_RULE
mkdir -p $(@D)
python $(SCRIPTDIR)/process-reads/cut-up-fasta.py $^ > $(@D)/velvet-noscaf-cut-up.fasta
runAssembly -force -o $(@D) $(@D)/velvet-noscaf-cut-up.fasta
rm $(@D)/velvet-noscaf-cut-up.fasta
mv $(@D)/454AllContigs.fna $@
endef
$(NEWBLER_OUT_VELVET_NOSCAF)/$(MERGE_FILENAME): $(VELVETG_OUT_NOSCAF)
	$(NEWBLER_RULE)
$(NEWBLER_OUT_METAVELVET_NOSCAF)/$(MERGE_FILENAME): $(METAVELVETG_OUT_NOSCAF)
	$(NEWBLER_RULE)
$(NEWBLER_OUT_RAY_NOSCAF)/$(MERGE_FILENAME): $(RAY_CONTIGS_OUT)
	$(NEWBLER_RULE)
################################
# -------- /newbler -----------#
################################

################################
# --------- bambus2 -----------#
################################
# Bambus2
define BAMBUS2_RULE
mkdir -p $(@D)
bash -x $(SCRIPTDIR)/map/map-bwa-markduplicates.sh $(MAP_PARS) $(FASTQ_TRIM_1) $(FASTQ_TRIM_2) \
	$(FASTQBASE) $< contigs $(@D)
bash -x $(SCRIPTDIR)/assembly/scaf-asm-bambus2.sh \
	$(@D)/contigs_${FASTQBASE}-smds.bam $< bambus2
endef
%/bambus2/bambus2.scaffold.linear.fasta: %/$(MERGE_FILENAME)
	$(BAMBUS2_RULE)
%/bambus2/bambus2.scaffold.linear.fasta: %/$(CONTIG_FILENAME)
	$(BAMBUS2_RULE)
################################
# -------- /bambus2 -----------#
################################

################################
# ----------- ray -------------#
################################
# Run Ray, do not use checkpoints, changes when nr of cores changes etc
define RAY_RULE
rm -rf $(@D)
$(MPI_EXEC_CMD) Ray -k $* -i $< -o $(@D) $(EXTRA_RAY_PARAMETERS)
endef
# Create symbolic link for pair, name must end on fastq for Ray
$(RAY_OUT)/pair.fastq: $(FASTQ_TRIM_IL)
	mkdir -p $(@D)
	ln -fs $(abspath $<) $@
$(RAY_OUT)/out_%/Contigs.fasta: $(RAY_OUT)/pair.fastq
	$(RAY_RULE)
$(RAY_OUT)/noscaf/noscaf_%/$(CONTIG_FILENAME): $(RAY_OUT)/out_%/Contigs.fasta
	mkdir -p $(@D)
	cp $(dir $<)Contigs.fasta $@
# Create links to the scaffold in a sub directory so validation is easier i.e. all output in a different folder
$(RAY_OUT)/scaf/scaf_%/$(SCAF_FILENAME): $(RAY_OUT)/noscaf/noscaf_%/$(CONTIG_FILENAME)
	mkdir -p $(@D)
	mv $(RAY_OUT)/out_$*/Scaffolds.fasta $@
	touch $@
################################
# ---------- /ray -------------#
################################
