################################
# ------ input parameters ---- #
################################
FASTQ1?=pair1.fastq
FASTQ2?=pair2.fastq
FASTQBASE?=pair
################################
# ----- /input parameters ---- #
################################

################################
# ----- output parameters ---- #
################################
# folder/file names
OUT?=ma-out
PRC_READS_OUT?=$(OUT)/processed-reads
ASM_OUT?=$(OUT)/assemblies
CONTIG_FILENAME?=ma-contigs.fa
SCAF_FILENAME?=ma-scaffolds.fa
MERGE_FILENAME?=ma-merge.fa
# K values
KMIN?=19
KMAX?=75
STEPSIZE?=2
KNUMBERS=$(shell seq $(KMIN) $(STEPSIZE) $(KMAX))
################################
# ---- /output parameters ---- #
################################

################################
# ------ quality trim -------- #
################################
DO_QTRIM?=yes
ifeq ($(DO_QTRIM),yes)
FASTQ_TRIM_1:=$(PRC_READS_OUT)/$(FASTQBASE).1.qtrim
FASTQ_TRIM_2:=$(PRC_READS_OUT)/$(FASTQBASE).2.qtrim
FASTQ_TRIM_IL:=$(PRC_READS_OUT)/$(FASTQBASE).qtrim
FASTQ_TRIM_UN:=$(FASTQ_TRIM_IL).unpaired
else
FASTQ_TRIM_1:=$(FASTQ1)
FASTQ_TRIM_2:=$(FASTQ2)
FASTQ_TRIM_IL?=$(PRC_READS_OUT)/$(FASTQBASE).qtrim
endif
FASTQ_TRIM_OUT:=$(FASTQ_TRIM_1) \
                $(FASTQ_TRIM_2) \
                $(FASTQ_TRIM_UN)
################################
# ------ /quality trim ------- #
################################

################################
# ---------- velveth --------- #
################################
VELVET_OUT:=$(ASM_OUT)/velvet
VELVETH_OUT:=$(VELVET_OUT)/velveth
VELVET_OUT_NOSCAF:=$(VELVET_OUT)/noscaf
VELVETH_OUT_SEQ:=$(foreach i,$(KNUMBERS),$(VELVETH_OUT)/velveth_$(i)/Sequences)
VELVETH_OUT_RD:=$(foreach i,$(KNUMBERS),$(VELVETH_OUT)/velveth_$(i)/Roadmaps)
VELVETG_OUT_NOSCAF:=$(foreach i,$(KNUMBERS),$(VELVET_OUT_NOSCAF)/noscaf_$(i)/$(CONTIG_FILENAME))
VELVET_OUT_SCAF:=$(VELVET_OUT)/scaf
VELVETG_OUT_SCAF:=$(foreach i,$(KNUMBERS),$(VELVET_OUT_SCAF)/scaf_$(i)/$(SCAF_FILENAME))
################################
# --------- /velveth --------- #
################################

################################
# ------- meta-velvetg ------- #
################################
METAVELVET_OUT:=$(ASM_OUT)/metavelvet
METAVELVET_OUT_NOSCAF:=$(METAVELVET_OUT)/noscaf
METAVELVETH_OUT_NOSCAF:=$(foreach i,$(KNUMBERS),$(METAVELVET_OUT_NOSCAF)/noscaf_$(i)/Sequences)
METAVELVETG_OUT_NOSCAF:=$(foreach i,$(KNUMBERS),$(METAVELVET_OUT_NOSCAF)/noscaf_$(i)/$(CONTIG_FILENAME))
METAVELVET_OUT_SCAF:=$(METAVELVET_OUT)/scaf
METAVELVETH_OUT_SCAF:=$(foreach i,$(KNUMBERS),$(METAVELVET_OUT_SCAF)/scaf_$(i)/Sequences)
METAVELVETG_OUT_SCAF:=$(foreach i,$(KNUMBERS),$(METAVELVET_OUT_SCAF)/scaf_$(i)/$(SCAF_FILENAME))
################################
# ------- /meta-velvetg ------ #
################################

################################
# ----------- ray -------------#
################################
MPI_EXEC_CMD?=mpiexec
EXTRA_RAY_PARAMETERS?=-show-memory-usage 
RAY_OUT:=$(ASM_OUT)/ray
RAY_OUT_NOSCAF:=$(RAY_OUT)/noscaf
RAY_OUT_SCAF:=$(RAY_OUT)/scaf
RAY_CONTIGS_OUT:=$(foreach i,$(KNUMBERS),$(RAY_OUT_NOSCAF)/noscaf_$(i)/$(CONTIG_FILENAME))
RAY_SCAFFOLDS_OUT:=$(foreach i,$(KNUMBERS),$(RAY_OUT_SCAF)/scaf_$(i)/$(SCAF_FILENAME))
################################
# ---------- /ray -------------#
################################

################################
# --------- minimus2  -------- #
################################
# Minimus2
MINIMUS2_OUT_VELVET_NOSCAF:=$(VELVET_OUT_NOSCAF)/minimus2
MINIMUS2_OUT_METAVELVET_NOSCAF:=$(METAVELVET_OUT_NOSCAF)/minimus2
MINIMUS2_OUT_RAY_NOSCAF:=$(RAY_OUT_NOSCAF)/minimus2
MINIMUS2_OUT_VELVET_MERGE:=$(MINIMUS2_OUT_VELVET_NOSCAF)/$(MERGE_FILENAME)
MINIMUS2_OUT_METAVELVET_MERGE:=$(MINIMUS2_OUT_METAVELVET_NOSCAF)/$(MERGE_FILENAME)
MINIMUS2_OUT_RAY_MERGE:=$(MINIMUS2_OUT_RAY_NOSCAF)/$(MERGE_FILENAME)
MINIMUS2_OUT_MERGE:=$(MINIMUS2_OUT_VELVET_MERGE) $(MINIMUS2_OUT_METAVELVET_MERGE) $(MINIMUS2_OUT_RAY_MERGE)
################################
# --------- /minimus2  ------- #
################################

################################
# --------- newbler -----------#
################################
NEWBLER_OUT_VELVET_NOSCAF:=$(VELVET_OUT_NOSCAF)/newbler
NEWBLER_OUT_METAVELVET_NOSCAF:=$(METAVELVET_OUT_NOSCAF)/newbler
NEWBLER_OUT_RAY_NOSCAF:=$(RAY_OUT_NOSCAF)/newbler
NEWBLER_OUT_VELVET_MERGE:=$(NEWBLER_OUT_VELVET_NOSCAF)/$(MERGE_FILENAME)
NEWBLER_OUT_METAVELVET_MERGE:=$(NEWBLER_OUT_METAVELVET_NOSCAF)/$(MERGE_FILENAME)
NEWBLER_OUT_RAY_MERGE:=$(NEWBLER_OUT_RAY_NOSCAF)/$(MERGE_FILENAME)
NEWBLER_OUT_MERGE:=$(NEWBLER_OUT_VELVET_MERGE) $(NEWBLER_OUT_METAVELVET_MERGE) $(NEWBLER_OUT_RAY_MERGE)
################################
# -------- /newbler -----------#
################################

################################
# --------- bambus2 -----------#
################################
BAMBUS2_MAP_PARS?=-t $(shell nproc)
################################
# -------- /bambus2 -----------#
################################

################################
# contig & scaffold filenames  #
################################
# A list of the full paths of all the contigs that can be assembled
ALLASMCONTIGS:=$(VELVETG_OUT_NOSCAF) \
			  $(METAVELVETG_OUT_NOSCAF) \
			  $(RAY_CONTIGS_OUT) \
			  $(MINIMUS2_OUT_MERGE) \
			  $(NEWBLER_OUT_MERGE) 
BAMBUS2SCAFFOLDS=$(foreach contigs, $(ALLASMCONTIGS), $(shell dirname $(contigs))/bambus2/bambus2.scaffold.linear.fasta)
BAMBUS2SCAFFOLDS_EXISTING=$(foreach contigs, $(wildcard $(ALLASMCONTIGS)), $(shell dirname $(contigs))/bambus2/bambus2.scaffold.linear.fasta)
ALLASMSCAFFOLDS=$(VELVETG_OUT_SCAF) $(METAVELVETG_OUT_SCAF) $(BAMBUS2SCAFFOLDS) $(RAY_SCAFFOLDS_OUT)
# Prints the full the paths of the assemblies files to stdout
echoall:
	@echo $(ALLASMCONTIGS) $(ALLASMSCAFFOLDS)
echocontigs:
	@echo $(ALLASMCONTIGS)
echoscaffolds:
	@echo $(ALLASMSCAFFOLDS)
echoexisting:
	@echo $(wildcard $(ALLASMCONTIGS) $(ALLASMSCAFFOLDS))
################################
# /contig & scaffold filenames #
################################

################################
#    Rules to make assemblies  #
################################
all: velvet metavelvet ray minimus2 newbler bambus2
qtrim: $(FASTQ_TRIM_IL)
velvet: $(VELVETG_OUT_NOSCAF) $(VELVETG_OUT_SCAF)
metavelvet: $(METAVELVETG_OUT_NOSCAF) $(METAVELVETG_OUT_SCAF)
minimus2: $(MINIMUS2_OUT_MERGE)
newbler: $(NEWBLER_OUT_MERGE)
velvetnoscafnewbler: $(NEWBLER_OUT_VELVET_MERGE)
ray: $(RAY_CONTIGS_OUT) $(RAY_SCAFFOLDS_OUT)
bambus2: $(BAMBUS2SCAFFOLDS)
bambus2existing: $(BAMBUS2SCAFFOLDS_EXISTING)
################################
#   /Rules to make assemblies  #
################################

################################
#  Rules to delete assemblies  #
################################
# Remove output, often one is only interested in keeping the assemblies, make keepcontigsonly allows you to do that
#keepcontigsonly:
#	-find $(ASM_OUT)/* -type f | grep -v $(foreach contigs,$(wildcard $(ALLASMCONTIGS)),-e $(contigs)) -e $(VELVET_OUT_NOSCAF)/noscaf_$(KMAX)/Sequences -e slurm | xargs rm
#clean:
#	-rm -rf $(OUT)
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
	-rm -rf $(MINIMUS2_OUT_VELVET_NOSCAF) $(MINIMUS2_OUT_METAVELVET_NOSCAF) $(MINIMUS2_OUT_RAY_NOSCAF)
cleannewbler:
	-rm -rf $(NEWBLER_OUT_METAVELVET_NOSCAF) $(NEWBLER_OUT_VELVET_NOSCAF) $(NEWBLER_OUT_RAY_NOSCAF)
################################
# /Rules to delete assemblies  #
################################

.PHONY: all qtrim velvet metavelvet cleanall cleanasm cleanvelvetg cleanvelvet cleanmetavelvet cleanmetavelvetg cleanqtrim cleanminimus2 cleannewbler validateexisting keepcontigsonly echoexisting ray bambus2existing bambus2
# Takes quite some time to compute some of these, so you might want to decide yourself when to delete them by using make keepcontigsonly for instance.
.PRECIOUS: $(VELVETH_OUT_RD) $(VELVETH_OUT_SEQ) $(FASTQ_TRIM_IL)
