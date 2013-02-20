# MetAssemble Pipeline. This Makefile is NOT meant to build the scripts but to
# run the metassemble scripts on an Illumina paired end library. Change input
# parameters accordingly for each library, either in the file or pass as
# argument. Currently only works with sbatch, that we use at our server. This
# version includes proper scheduling of job dependencies except for bambus2. I
# do those separately with 'make bambus2' after constructing the contigs
# because I have not found time to implement that yet and a lot of assemblies
# fail anyway so it didn't make much sense to schedule the bambus2 jobs
# immediately. The Makefile can also be run without SBATCH by setting
# USE_SBATCH=no. Please use 'make -j NROFJOBS' in that case to parallelize. If
# one wants the commands for a specific assembly one can use 'make -Bn
# path/to/contigs.fasta', e.g.:
# 'make -Bn ma-out/velvet/noscaf/noscaf_31/contigs.fasta'
include /gulo/glob/inod/github/metassemble/scripts/parameters.mk

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
# qual-type 0 for Sanger quality (or illumina 1.8) type 2 for windowed trimming
$(PRC_READS_OUT)/$(FASTQBASE).qtrim: $(FASTQ1) $(FASTQ2)
	mkdir -p $(PRC_READS_OUT)
	sickle pe \
		-f $(FASTQ1) \
		-r $(FASTQ2) \
		-t sanger \
		-o $(PRC_READS_OUT)/$(FASTQBASE).1.qtrim \
		-p $(PRC_READS_OUT)/$(FASTQBASE).2.qtrim \
		-s $(PRC_READS_OUT)/$(FASTQBASE).qtrim.unpaired
	shuffleSequences_fastq.pl $(PRC_READS_OUT)/$(FASTQBASE).1.qtrim \
		$(PRC_READS_OUT)/$(FASTQBASE).2.qtrim \
		$(PRC_READS_OUT)/$(FASTQBASE).qtrim
################################
# ------ /quality trim ------- #
################################

################################
# ---------- velveth --------- #
################################
$(VELVETH_OUT)/velveth_$(KMIN)/Sequences: $(PRC_READS_OUT)/$(FASTQBASE).qtrim
	mkdir -p $(VELVETH_OUT)
	velveth $(VELVETH_OUT)/velveth_$(KMIN) $(KMIN) -noHash -fastq \
		-shortPaired $(PRC_READS_OUT)/$(FASTQBASE).qtrim
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
endef
$(VELVET_OUT_NOSCAF)/noscaf_%/contigs.fa: $(VELVETH_OUT)/velveth_%/Roadmaps $(VELVETH_OUT)/velveth_%/Sequences
	$(call velvetg_rule,-scaffolding no)
$(VELVET_OUT_SCAF)/scaf_%/contigs.fa: $(VELVETH_OUT)/velveth_%/Roadmaps $(VELVETH_OUT)/velveth_%/Sequences
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
ln -fs $(lastword $^) $(@D)/Sequences
ln -fs $< $(@D)/Roadmaps
velvetg $(dir $@) -scaffolding no -exp_cov auto -read_trkg yes \
	&& meta-velvetg $(dir $@) $1
endef
$(METAVELVET_OUT_NOSCAF)/noscaf_%/meta-velvetg.contigs.fa: $(VELVETH_OUT)/velveth_%/Roadmaps $(VELVETH_OUT)/velveth_%/Sequences
	$(call metavelvetg_rule,-scaffolding no)
$(METAVELVET_OUT_SCAF)/scaf_%/meta-velvetg.contigs.fa: $(VELVETH_OUT)/velveth_%/Roadmaps $(VELVETH_OUT)/velveth_%/Sequences
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
bash $(SCRIPTDIR)/assembly/merge-asm-minimus2.sh $(@D) $^
endef
$(MINIMUS2_OUT_VELVET_NOSCAF)/all-merged.fasta: $(VELVETG_OUT_NOSCAF)
	$(MINIMUS2_RULE)
$(MINIMUS2_OUT_METAVELVET_NOSCAF)/all-merged.fasta: $(METAVELVETG_OUT_NOSCAF)
	$(MINIMUS2_RULE)
$(MINIMUS2_OUT_RAY_NOSCAF)/all-merged.fasta: $(RAY_CONTIGS_OUT)
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
endef
$(NEWBLER_OUT_VELVET_NOSCAF)/454AllContigs.fna: $(VELVETG_OUT_NOSCAF)
	$(NEWBLER_RULE)
$(NEWBLER_OUT_METAVELVET_NOSCAF)/454AllContigs.fna: $(METAVELVETG_OUT_NOSCAF)
	$(NEWBLER_RULE)
$(NEWBLER_OUT_RAY_NOSCAF)/454AllContigs.fna: $(RAY_CONTIGS_OUT)
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
bash $(SCRIPTDIR)/validate/map-bwa-markduplicates.sh $(FASTQ1) $(FASTQ2) \
	$(FASTQBASE) $< contigs $(@D)
bash $(SCRIPTDIR)/assembly/scaf-asm-bambus2.sh \
	$(@D)/contigs_${FASTQBASE}-smds.bam $< bambus2
rm $(@D)/contigs_${FASTQBASE}-smds.bam
endef
%/bambus2/bambus2.scaffold.linear.fasta: %/all-merged.fasta
	$(BAMBUS2_RULE)
%/bambus2/bambus2.scaffold.linear.fasta: %/contigs.fa
	$(BAMBUS2_RULE)
%/bambus2/bambus2.scaffold.linear.fasta: %/454AllContigs.fna
	$(BAMBUS2_RULE)
%/bambus2/bambus2.scaffold.linear.fasta: %/meta-velvetg.contigs.fa
	$(BAMBUS2_RULE)
%/bambus2/bambus2.scaffold.linear.fasta: %/ma-contigs.fasta
	$(BAMBUS2_RULE)
################################
# -------- /bambus2 -----------#
################################

################################
# ----------- ray -------------#
################################
define RAY_RULE
rm -rf $(@D) $(@D).cp
mpiexec Ray -k $* -i $< -o $(@D) -show-memory-usage -write-checkpoints $(@D).cp
endef
# Create symbolic link for pair, name must end on fastq for Ray
$(RAY_OUT)/pair.fastq: $(PRC_READS_OUT)/$(FASTQBASE).qtrim
	mkdir -p $(@D)
	ln -fs $(abspath $<) $@
$(RAY_OUT)/out_%/Contigs.fasta: $(RAY_OUT)/pair.fastq
	$(RAY_RULE)
$(RAY_OUT)/noscaf/noscaf_%/ma-contigs.fasta: $(RAY_OUT)/out_%/Contigs.fasta
	mkdir -p $(@D)
	mv $(dir $<)Contigs.fasta $@
# Create links to the scaffold in a sub directory so validation is easier i.e. all output in a different folder
$(RAY_OUT)/scaf/scaf_%/ma-scaffolds.fasta: $(RAY_OUT)/out_%/Contigs.fasta
	mkdir -p $(@D)
	mv $(dir $<)Scaffolds.fasta $@
################################
# ---------- /ray -------------#
################################
