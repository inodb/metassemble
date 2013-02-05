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

################################
# ------ input parameters ---- #
################################
FASTQ1=500pg_even_GGACTCCT-CTCTCTAT_L003_R1_001.fastq
FASTQ2=500pg_even_GGACTCCT-CTCTCTAT_L003_R2_001.fastq
FASTQBASE=500pg_even
SCRIPTDIR=/bubo/home/h16/inod/glob/github/metassemble/scripts
# sbatch specific input
USE_SBATCH=no
SB_ACCOUNT=b2010008
SB_USER=inod
SB_WRAPPER_SCRIPT=/bubo/home/h16/inod/bin/sbatch_job # like --wrap with extras
################################
# ----- /input parameters ---- #
################################

################################
# ----- output parameters ---- #
################################
OUT=ma-out
PRC_READS_OUT=$(OUT)/processed-reads
ASM_OUT=$(OUT)/assemblies
################################
# ---- /output parameters ---- #
################################


################################
# ------ sbatch functions ---- #
################################
# Function for some default parameters for sbatch job, the shell command
# returns nothing, it is just to make sure output directory exists for the
# slurm file.
# $1 is rulename, $3 is directory and $2 are sbatch parameters
get_sbatch_job_par=-J $1.$(FASTQBASE) --output=$3slurm-$1.out --error=$3slurm-$1.err \
				   -A $(SB_ACCOUNT) $2 $(SB_WRAPPER_SCRIPT) \
				   $(shell mkdir -p $3) 
# Generates the job dependency parameters - works only for two job names atm,
# but can be easily extended, not sure how to generalize for all given
# arguments
define get_job_dependencies
$(eval jobnames=$1.$(FASTQBASE))
$(if $2,$(eval jobnames=$(jobnames),$2.$(FASTQBASE)))
squeue --noheader -u $(SB_USER) -n $(jobnames) | awk '{if (NR == 1) {printf "-d afterok"} printf(":%s",$$1)}'
endef
################################
# ----- /sbatch functions ---- #
################################

################################
# ------ quality trim -------- #
################################
FASTQ_TRIM_OUT=$(PRC_READS_OUT)/$(FASTQBASE).1.qtrim \
               $(PRC_READS_OUT)/$(FASTQBASE).2.qtrim \
               $(PRC_READS_OUT)/$(FASTQBASE).qtrim.unpaired
# qual-type 0 for Sanger quality (or illumina 1.8) type 2 for windowed trimming
QTRIM_RULE=\
    mkdir -p $(PRC_READS_OUT); \
    sickle pe \
        -f $(FASTQ1) \
        -r $(FASTQ2) \
        -t sanger \
        -o $(PRC_READS_OUT)/$(FASTQBASE).1.qtrim \
        -p $(PRC_READS_OUT)/$(FASTQBASE).2.qtrim \
        -s $(PRC_READS_OUT)/$(FASTQBASE).qtrim.unpaired; \
    shuffleSequences_fastq.pl $(PRC_READS_OUT)/$(FASTQBASE).1.qtrim \
        $(PRC_READS_OUT)/$(FASTQBASE).2.qtrim \
        $(PRC_READS_OUT)/$(FASTQBASE).qtrim
$(PRC_READS_OUT)/$(FASTQBASE).qtrim : $(FASTQ1) $(FASTQ2)
ifeq ($(USE_SBATCH),yes)
	sbatch $(call get_sbatch_job_par,qtrim,-p core -t 6:00:00,$(@D)) "$(QTRIM_RULE)"
else
	$(QTRIM_RULE)
endif
################################
# ------ /quality trim ------- #
################################

################################
# ---------- velveth --------- #
################################
KMIN=31
KMAX=84
STEPSIZE=2
# Output OKMAX < KMAX in velveth, so it is actually like
# {KMIN..(KMAX-1)..STEPSIZE} in bash
OKMAX=$(shell echo $$(($(KMAX) - 1 - (($(KMAX) - 1 - $(KMIN)) % $(STEPSIZE)))))
VELVET_OUT=$(ASM_OUT)/velvet
VELVETH_OUT=$(VELVET_OUT)/velveth
VELVET_OUT_NOSCAF=$(VELVET_OUT)/noscaf
VELVETH_OUT_SIGNAL_FILE=$(VELVETH_OUT)/velveth.done
VELVETH_OUT_SEQ=$(shell echo $(VELVETH_OUT)/velveth_{$(KMIN)..$(OKMAX)..$(STEPSIZE)}/Sequences)
VELVETG_OUT_NOSCAF=$(shell echo $(VELVET_OUT_NOSCAF)/noscaf_{$(KMIN)..$(OKMAX)..$(STEPSIZE)}/contigs.fa)
VELVET_OUT_SCAF=$(VELVET_OUT)/scaf
VELVETG_OUT_SCAF=$(shell echo $(VELVET_OUT_SCAF)/scaf_{$(KMIN)..$(OKMAX)..$(STEPSIZE)}/contigs.fa)
# velveth
VELVETH_RULE=\
    mkdir -p $(VELVETH_OUT); \
    velveth $(VELVETH_OUT)/velveth $(KMIN),$(KMAX),$(STEPSIZE) \
        -fastq -shortPaired $(PRC_READS_OUT)/$(FASTQBASE).qtrim; \
	touch $(VELVETH_OUT_SIGNAL_FILE)
$(VELVETH_OUT_SIGNAL_FILE) : $(PRC_READS_OUT)/$(FASTQBASE).qtrim
ifeq ($(USE_SBATCH),yes)
	sbatch $(shell $(call get_job_dependencies,qtrim)) $(call get_sbatch_job_par,velveth,-p node -C mem24GB -t 1-00:00:00,$(@D)) \
		"$(VELVETH_RULE)"
else
	$(VELVETH_RULE)
endif
$(VELVETH_OUT)/velveth_%/Sequences : $(VELVETH_OUT_SIGNAL_FILE)
################################
# --------- /velveth --------- #
################################

################################
# --------- velvetg ---------- #
################################
velvetg_rule=\
    mkdir -p $(dir $@); \
    velvetg $(dir $<) $1 || { echo 'velvetg $1 failed.' >&2; }; \
    mv $(dir $<)contigs.fa $(dir $@)
# No scaffolding velvetg, perform on velveth output and store output in a
# different folder, so velveth output can be reused
$(VELVET_OUT_NOSCAF)/noscaf_%/contigs.fa : $(VELVETH_OUT)/velveth_%/Sequences
ifeq ($(USE_SBATCH),yes)
	sbatch $(shell $(call get_job_dependencies,velveth)) $(call get_sbatch_job_par,vgnoscaf$*,-p node -t 02:00:00,$(@D)) \
		"$(call velvetg_rule,-scaffolding no)"
else
	$(call velvetg_rule,-scaffolding no)
endif
# Scaffolding velvetg, -exp_cov auto although unsuitable for metagenomics must
# be set otherwise scaffolding is not performed
$(VELVET_OUT_SCAF)/scaf_%/contigs.fa : $(VELVETH_OUT)/velveth_%/Sequences
ifeq ($(USE_SBATCH),yes)
	sbatch $(shell $(call get_job_dependencies,velveth)) $(call get_sbatch_job_par,vgscaf$*,-p node -t 02:00:00,$(@D)) \
		"$(call velvetg_rule,-scaffolding no -exp_cov auto)"
else
	$(call velvetg_rule,-scaffolding no -exp_cov auto)
endif
################################
# --------- /velvetg --------- #
################################

################################
# ------- meta-velvetg ------- #
################################
METAVELVET_OUT=$(ASM_OUT)/metavelvet
METAVELVET_OUT_NOSCAF=$(METAVELVET_OUT)/noscaf
METAVELVETH_OUT_NOSCAF=$(shell echo $(METAVELVET_OUT_NOSCAF)/noscaf_{$(KMIN)..$(OKMAX)..$(STEPSIZE)}/Sequences)
METAVELVETG_OUT_NOSCAF=$(shell echo $(METAVELVET_OUT_NOSCAF)/noscaf_{$(KMIN)..$(OKMAX)..$(STEPSIZE)}/meta-velvetg.contigs.fa)
METAVELVET_OUT_SCAF=$(METAVELVET_OUT)/scaf
METAVELVETH_OUT_SCAF=$(shell echo $(METAVELVET_OUT_SCAF)/scaf_{$(KMIN)..$(OKMAX)..$(STEPSIZE)}/Sequences)
METAVELVETG_OUT_SCAF=$(shell echo $(METAVELVET_OUT_SCAF)/scaf_{$(KMIN)..$(OKMAX)..$(STEPSIZE)}/meta-velvetg.contigs.fa)
# Copy output from velveth and run velvetg, followed by meta-velvetg -scaffolding yes or no
metavelvetg_rule=\
    mkdir -p $(dir $@); \
    velvetg $(dir $<) -scaffolding no -exp_cov auto -read_trkg yes \
    && meta-velvetg $(dir $<) $1 \
    || { echo 'velvetg or meta-velvetg $1 failed.' >&2; }; \
    mv $(dir $<)meta-velvetg.contigs.fa $(dir $@)
$(METAVELVET_OUT_NOSCAF)/noscaf_%/meta-velvetg.contigs.fa : $(VELVETH_OUT)/velveth_%/Sequences
ifeq ($(USE_SBATCH),yes)
	sbatch $(shell $(call get_job_dependencies,velveth)) $(call get_sbatch_job_par,mvgnoscaf$*,-p node -t 02:00:00,$(@D)) \
		"$(call metavelvetg_rule,-scaffolding no)"
else
	$(call metavelvetg_rule,-scaffolding no)
endif
$(METAVELVET_OUT_SCAF)/scaf_%/meta-velvetg.contigs.fa : $(VELVETH_OUT)/velveth_%/Sequences
ifeq ($(USE_SBATCH),yes)
	sbatch $(shell $(call get_job_dependencies,velveth)) $(call get_sbatch_job_par,mvgscaf$*,-p node -t 02:00:00,$(@D)) \
		"$(call metavelvetg_rule,-scaffolding yes)"
else
	$(call metavelvetg_rule,-scaffolding yes)
endif
################################
# ------- /meta-velvetg ------ #
################################

################################
# --------- minimus2  -------- #
################################
# Minimus2
MINIMUS2_OUT_VELVET_NOSCAF=$(VELVET_OUT_NOSCAF)/minimus2
MINIMUS2_OUT_METAVELVET_NOSCAF=$(METAVELVET_OUT_NOSCAF)/minimus2
# Minimus2 rule merges all given prerequisite files
MINIMUS2_RULE=\
    mkdir -p $(dir $@); \
    bash $(SCRIPTDIR)/assembly/merge-asm-minimus2.sh \
        $(@D) $^
$(MINIMUS2_OUT_VELVET_NOSCAF)/all-merged.fasta : $(VELVETG_OUT_NOSCAF)
ifeq ($(USE_SBATCH),yes)
	sbatch $(shell $(call get_job_dependencies,vgnoscaf$(OKMAX))) $(call get_sbatch_job_par,minimus2vnoscaf,-p node -t 1-00:00:00,$(@D)) \
		"$(MINIMUS2_RULE)"
else
	$(MINIMUS2_RULE)
endif
$(MINIMUS2_OUT_METAVELVET_NOSCAF)/all-merged.fasta : $(METAVELVETG_OUT_NOSCAF)
ifeq ($(USE_SBATCH),yes)
	sbatch $(shell $(call get_job_dependencies,mvgnoscaf$(OKMAX))) $(call get_sbatch_job_par,minimus2mvnoscaf,-p node -t 1-00:00:00,$(@D)) \
		"$(MINIMUS2_RULE)"
else
	$(MINIMUS2_RULE)
endif
################################
# --------- /minimus2  ------- #
################################

################################
# --------- newbler -----------#
################################
NEWBLER_OUT_VELVET_NOSCAF=$(VELVET_OUT_NOSCAF)/newbler
NEWBLER_OUT_METAVELVET_NOSCAF=$(METAVELVET_OUT_NOSCAF)/newbler
# Newbler rule merges all given prerequisite files
NEWBLER_RULE=\
    mkdir -p $(dir $@); \
    module load biopython; \
    python $(SCRIPTDIR)/process-reads/cut-up-fasta.py \
        $^ > $(@D)/velvet-noscaf-cut-up.fasta; \
    module load 454-dataanalysis/2.6; \
    runAssembly -force -o $(@D) \
        $(@D)/velvet-noscaf-cut-up.fasta; \
    rm $(@D)/velvet-noscaf-cut-up.fasta
$(NEWBLER_OUT_VELVET_NOSCAF)/454AllContigs.fna: $(VELVETG_OUT_NOSCAF)
ifeq ($(USE_SBATCH),yes)
	sbatch $(shell $(call get_job_dependencies,vgnoscaf$(OKMAX))) $(call get_sbatch_job_par,newblervnoscaf,-p node -t 1-00:00:00,$(@D)) \
		"$(NEWBLER_RULE)"
else
	$(NEWBLER_RULE)
endif
$(NEWBLER_OUT_METAVELVET_NOSCAF)/454AllContigs.fna: $(METAVELVETG_OUT_NOSCAF)
ifeq ($(USE_SBATCH),yes)
	sbatch $(shell $(call get_job_dependencies,mvgnoscaf$(OKMAX))) $(call get_sbatch_job_par,newblermvnoscaf,-p node -t 1-00:00:00,$(@D)) \
		"$(NEWBLER_RULE)"
else
	$(NEWBLER_RULE)
endif
################################
# -------- /newbler -----------#
################################

################################
# --------- bambus2 -----------#
################################
# Bambus2
BAMBUS2_RULE=\
	mkdir -p $(@D); \
	bash $(SCRIPTDIR)/validate/map-bwa-markduplicates.sh \
		$(FASTQ1) $(FASTQ2) $(FASTQBASE) $< contigs $(@D); \
	bash $(SCRIPTDIR)/assembly/scaf-asm-bambus2.sh \
		$(@D)/contigs_${FASTQBASE}-smds.bam \
		$< bambus2; \
	rm $(@D)/contigs_${FASTQBASE}-smds.bam; 
%/bambus2/bambus2.scaffold.linear.fasta: %/all-merged.fasta
ifeq ($(USE_SBATCH),yes)
	sbatch $(call get_sbatch_job_par,bambus2minimus2,-p node -t 1-00:00:00,$(@D)) \
	"$(BAMBUS2_RULE)"
else
	$(BAMBUS2_RULE)
endif
%/bambus2/bambus2.scaffold.linear.fasta: %/contigs.fa
ifeq ($(USE_SBATCH),yes)
	sbatch $(call get_sbatch_job_par,bambus2vnoscaf,-p node -t 1-00:00:00,$(@D)) \
	"$(BAMBUS2_RULE)"
else
	$(BAMBUS2_RULE)
endif
%/bambus2/bambus2.scaffold.linear.fasta: %/454AllContigs.fna
ifeq ($(USE_SBATCH),yes)
	sbatch $(call get_sbatch_job_par,bambus2newbler,-p node -t 1-00:00:00,$(@D)) \
	"$(BAMBUS2_RULE)"
else
	$(BAMBUS2_RULE)
endif
%/bambus2/bambus2.scaffold.linear.fasta: %/meta-velvetg.contigs.fa
ifeq ($(USE_SBATCH),yes)
	sbatch $(call get_sbatch_job_par,bambus2mvnoscaf,-p node -t 1-00:00:00,$(@D)) \
	"$(BAMBUS2_RULE)"
else
	$(BAMBUS2_RULE)
endif
################################
# -------- /bambus2 -----------#
################################

################################
# ----------- ray -------------#
################################
RAY_OUT=$(ASM_OUT)/ray
RAY_CONTIGS_OUT=$(shell echo $(RAY_OUT)/out_{$(KMIN)..$(KMAX)..$(STEPSIZE)}/Contigs.fasta)
RAY_SCAFFOLDS_OUT=$(shell echo $(RAY_OUT)/out_{$(KMIN)..$(KMAX)..$(STEPSIZE)}/scaf/Scaffolds.fasta)
# Loop iterates over single k
ray_rule=\
    rm -rf $(@D) $(@D).cp; \
	module load openmpi; \
	mpiexec Ray -k $* -i $< -o $(@D) -show-memory-usage -write-checkpoints $(@D).cp
# Create symbolic link for pair, name must end on fastq for Ray
$(RAY_OUT)/pair.fastq : $(PRC_READS_OUT)/$(FASTQBASE).qtrim
	ln -fs $(shell readlink -f $<) $@
$(RAY_OUT)/out_%/Contigs.fasta : $(RAY_OUT)/pair.fastq
ifeq ($(USE_SBATCH),yes)
	sbatch $(shell $(call get_job_dependencies,qtrim)) $(call get_sbatch_job_par,ray$*,-p node -t 12:00:00,$(RAY_OUT)/) \
		"$(call ray_rule)"
else
	$(call ray_rule)
endif
# Create links to the scaffold in a sub directory so validation is easier i.e. all output in a different folder
$(RAY_OUT)/out_%/scaf/Scaffolds.fasta : $(RAY_OUT)/out_%/Scaffolds.fasta
	mkdir -p $(@D)
	ln -fs $(shell readlink -f $<) $@
################################
# ---------- /ray -------------#
################################
	

################################
# contig & scaffold filenames  #
################################
# A list of the full paths of all the contigs that can be assembled
ALLASMCONTIGS=$(VELVETG_OUT_NOSCAF) \
			  $(METAVELVETG_OUT_NOSCAF) \
			  $(MINIMUS2_OUT_VELVET_NOSCAF)/all-merged.fasta \
			  $(MINIMUS2_OUT_METAVELVET_NOSCAF)/all-merged.fasta \
			  $(NEWBLER_OUT_VELVET_NOSCAF)/454AllContigs.fna \
			  $(NEWBLER_OUT_METAVELVET_NOSCAF)/454AllContigs.fna \
			  $(RAY_CONTIGS_OUT)
BAMBUS2SCAFFOLDS=$(foreach contigs, $(ALLASMCONTIGS), $(dir $(contigs))bambus2/bambus2.scaffold.linear.fasta)
ALLASMSCAFFOLDS=$(VELVETG_OUT_SCAF) $(METAVELVETG_OUT_SCAF) $(BAMBUS2SCAFFOLDS) $(RAY_SCAFFOLDS_OUT)
# Prints the full the paths of the assemblies files to stdout, useful e.g. if one wants to only copy
# contigs someplace else
echoexisting:
	@echo $(wildcard $(ALLASMCONTIGS) $(ALLASMSCAFFOLDS))
################################
# /contig & scaffold filenames #
################################

################################
#    Rules to make assemblies  #
################################
all: qtrim velvet metavelvet minimus2 newbler ray
qtrim: $(PRC_READS_OUT)/$(FASTQBASE).qtrim
velvet: $(VELVETG_OUT_NOSCAF) $(VELVETG_OUT_SCAF)
metavelvet: $(METAVELVETG_OUT_NOSCAF) $(METAVELVETG_OUT_SCAF)
ray: $(RAY_CONTIGS_OUT)
# The ray scaffolds are created automatically with ray, rayscaf just creates a
# symbolic link in a subdir so validation is easier (all output in its own dir)
rayscaf: $(RAY_SCAFFOLDS_OUT)

minimus2: minimus2velvet minimus2metavelvet
minimus2velvet: $(MINIMUS2_OUT_VELVET_NOSCAF)/all-merged.fasta
minimus2metavelvet: $(MINIMUS2_OUT_METAVELVET_NOSCAF)/all-merged.fasta

newbler: newblervelvet newblermetavelvet
newblervelvet: $(NEWBLER_OUT_VELVET_NOSCAF)/454AllContigs.fna
newblermetavelvet: $(NEWBLER_OUT_METAVELVET_NOSCAF)/454AllContigs.fna

bambus2: $(BAMBUS2SCAFFOLDS)

################################
#   /Rules to make assemblies  #
################################

################################
#  Rules to delete assemblies  #
################################
# Remove output, often one is only interested in keeping the assemblies, make keepcontigsonly allows you to do that
keepcontigsonly:
	-find $(ASM_OUT)/* -type f | grep -v $(foreach contigs,$(wildcard $(ALLASMCONTIGS)),-e $(contigs)) -e $(VELVET_OUT_NOSCAF)/noscaf_$(OKMAX)/Sequences -e slurm | xargs rm
clean:
	-rm -rf $(OUT)
cleanasm:
	-rm -rf $(ASM_OUT)
cleanvelvetg:
	-rm $(VELVETG_OUT_NOSCAF) $(VELVETG_OUT_SCAF)
cleanvelvet:
	-rm -rf $(VELVET_OUT)
cleanmetavelvet:
	-rm -rf $(METAVELVET_OUT)
cleanmetavelvetg:
	-rm $(METAVELVETG_OUT_NOSCAF) $(METAVELVETG_OUT_SCAF)
cleanqtrim:
	-rm -rf $(PRC_READS_OUT)
cleanminimus2:
	-rm -rf $(MINIMUS2_OUT_VELVET_NOSCAF) $(MINIMUS2_OUT_METAVELVET_NOSCAF)
cleannewbler:
	-rm -rf $(NEWBLER_OUT_METAVELVET_NOSCAF) $(NEWBLER_OUT_VELVET_NOSCAF)
cleanpurity:
	-find $(ASM_OUT)/* -type f -name purity-length-hist.tab | xargs rm
################################
# /Rules to delete assemblies  #
################################

.PHONY: all qtrim velvet metavelvet minimus2 minimus2velvet minimus2metavelvet newbler newblervelvet newblermetavelvet cleanall cleanasm cleanvelvetg cleanvelvet cleanmetavelvet cleanmetavelvetg cleanqtrim cleanminimus2 cleannewbler validateexisting keepcontigsonly echoexisting ray
# Takes quite some time to compute some of these, so you might want to decide yourself when to delete them by using make keepcontigsonly for instance.
.PRECIOUS: $(VELVETH_OUT_SEQ) $(PRC_READS_OUT)/$(FASTQBASE).qtrim
