################################
# ------ input parameters ---- #
################################
FASTQ1:=500pg_even_GGACTCCT-CTCTCTAT_L003_R1_001.fastq
FASTQ2:=500pg_even_GGACTCCT-CTCTCTAT_L003_R2_001.fastq
FASTQBASE:=500pg_even
SCRIPTDIR:=/bubo/home/h16/inod/glob/github/metassemble/scripts
################################
# ----- /input parameters ---- #
################################

################################
# ----- output parameters ---- #
################################
OUT:=ma-out
PRC_READS_OUT:=$(OUT)/processed-reads
ASM_OUT:=$(OUT)/assemblies
################################
# ---- /output parameters ---- #
################################

################################
# ------ quality trim -------- #
################################
FASTQ_TRIM_OUT:=$(PRC_READS_OUT)/$(FASTQBASE).1.qtrim \
                $(PRC_READS_OUT)/$(FASTQBASE).2.qtrim \
                $(PRC_READS_OUT)/$(FASTQBASE).qtrim.unpaired
################################
# ------ /quality trim ------- #
################################

################################
# ---------- velveth --------- #
################################
KMIN:=11
KMAX:=84
STEPSIZE:=2
VELVET_OUT:=$(ASM_OUT)/velvet
VELVETH_OUT:=$(VELVET_OUT)/velveth
VELVET_OUT_NOSCAF:=$(VELVET_OUT)/noscaf
VELVETH_OUT_SEQ:=$(shell echo $(VELVETH_OUT)/velveth_{$(KMIN)..$(KMAX)..$(STEPSIZE)}/Sequences)
VELVETH_OUT_RD:=$(shell echo $(VELVETH_OUT)/velveth_{$(KMIN)..$(KMAX)..$(STEPSIZE)}/Roadmaps)
VELVETG_OUT_NOSCAF:=$(shell echo $(VELVET_OUT_NOSCAF)/noscaf_{$(KMIN)..$(KMAX)..$(STEPSIZE)}/contigs.fa)
VELVET_OUT_SCAF:=$(VELVET_OUT)/scaf
VELVETG_OUT_SCAF:=$(shell echo $(VELVET_OUT_SCAF)/scaf_{$(KMIN)..$(KMAX)..$(STEPSIZE)}/contigs.fa)
################################
# --------- /velveth --------- #
################################

################################
# ------- meta-velvetg ------- #
################################
METAVELVET_OUT:=$(ASM_OUT)/metavelvet
METAVELVET_OUT_NOSCAF:=$(METAVELVET_OUT)/noscaf
METAVELVETH_OUT_NOSCAF:=$(shell echo $(METAVELVET_OUT_NOSCAF)/noscaf_{$(KMIN)..$(KMAX)..$(STEPSIZE)}/Sequences)
METAVELVETG_OUT_NOSCAF:=$(shell echo $(METAVELVET_OUT_NOSCAF)/noscaf_{$(KMIN)..$(KMAX)..$(STEPSIZE)}/meta-velvetg.contigs.fa)
METAVELVET_OUT_SCAF:=$(METAVELVET_OUT)/scaf
METAVELVETH_OUT_SCAF:=$(shell echo $(METAVELVET_OUT_SCAF)/scaf_{$(KMIN)..$(KMAX)..$(STEPSIZE)}/Sequences)
METAVELVETG_OUT_SCAF:=$(shell echo $(METAVELVET_OUT_SCAF)/scaf_{$(KMIN)..$(KMAX)..$(STEPSIZE)}/meta-velvetg.contigs.fa)
################################
# ------- /meta-velvetg ------ #
################################

################################
# ----------- ray -------------#
################################
RAY_OUT:=$(ASM_OUT)/ray
RAY_OUT_NOSCAF:=$(RAY_OUT)/noscaf
RAY_OUT_SCAF:=$(RAY_OUT)/scaf
RAY_CONTIGS_OUT:=$(shell echo $(RAY_OUT_NOSCAF)/noscaf_{$(KMIN)..$(KMAX)..$(STEPSIZE)}/ma-contigs.fasta)
RAY_SCAFFOLDS_OUT:=$(shell echo $(RAY_OUT_SCAF)/scaf_{$(KMIN)..$(KMAX)..$(STEPSIZE)}/ma-scaffolds.fasta)
################################
# ---------- /ray -------------#
################################

################################
# --------- minimus2  -------- #
################################
# Minimus2
MINIMUS2_OUT_VELVET_NOSCAF=$(VELVET_OUT_NOSCAF)/minimus2
MINIMUS2_OUT_METAVELVET_NOSCAF=$(METAVELVET_OUT_NOSCAF)/minimus2
MINIMUS2_OUT_RAY_NOSCAF=$(RAY_OUT_NOSCAF)/minimus2
################################
# --------- /minimus2  ------- #
################################

################################
# --------- newbler -----------#
################################
NEWBLER_OUT_VELVET_NOSCAF:=$(VELVET_OUT_NOSCAF)/newbler
NEWBLER_OUT_METAVELVET_NOSCAF:=$(METAVELVET_OUT_NOSCAF)/newbler
NEWBLER_OUT_RAY_NOSCAF:=$(RAY_OUT_NOSCAF)/newbler
################################
# -------- /newbler -----------#
################################

################################
# contig & scaffold filenames  #
################################
# A list of the full paths of all the contigs that can be assembled
ALLASMCONTIGS=$(VELVETG_OUT_NOSCAF) \
			  $(METAVELVETG_OUT_NOSCAF) \
			  $(MINIMUS2_OUT_VELVET_NOSCAF)/all-merged.fasta \
			  $(MINIMUS2_OUT_METAVELVET_NOSCAF)/all-merged.fasta \
			  $(MINIMUS2_OUT_RAY_NOSCAF)/all-merged.fasta \
			  $(NEWBLER_OUT_VELVET_NOSCAF)/454AllContigs.fna \
			  $(NEWBLER_OUT_METAVELVET_NOSCAF)/454AllContigs.fna \
			  $(NEWBLER_OUT_RAY_NOSCAF)/all-merged.fasta \
			  $(RAY_CONTIGS_OUT)
BAMBUS2SCAFFOLDS=$(foreach contigs, $(ALLASMCONTIGS), $(shell dirname $(contigs))/bambus2/bambus2.scaffold.linear.fasta)
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
all: qtrim velvet metavelvet ray minimus2 newbler
qtrim: $(PRC_READS_OUT)/$(FASTQBASE).qtrim
velvet: $(VELVETG_OUT_NOSCAF) $(VELVETG_OUT_SCAF)
metavelvet: $(METAVELVETG_OUT_NOSCAF) $(METAVELVETG_OUT_SCAF)
# Ray scaffolds are automatically made so contigs will be made as well
# NOTE: contigs has been removed for sbatch rule, otherwise you'll have to
# launch a seperate job to only create some symlinks.
ray: $(RAY_SCAFFOLDS_OUT)

minimus2: minimus2velvet minimus2metavelvet
minimus2velvet: $(MINIMUS2_OUT_VELVET_NOSCAF)/all-merged.fasta
minimus2metavelvet: $(MINIMUS2_OUT_METAVELVET_NOSCAF)/all-merged.fasta
minimus2ray: $(MINIMUS2_OUT_RAY_NOSCAF)/all-merged.fasta

newbler: newblervelvet newblermetavelvet
newblervelvet: $(NEWBLER_OUT_VELVET_NOSCAF)/454AllContigs.fna
newblermetavelvet: $(NEWBLER_OUT_METAVELVET_NOSCAF)/454AllContigs.fna
newblerray: $(NEWBLER_OUT_RAY_NOSCAF)/454AllContigs.fna

bambus2: $(BAMBUS2SCAFFOLDS)
################################
#   /Rules to make assemblies  #
################################

################################
#  Rules to delete assemblies  #
################################
# Remove output, often one is only interested in keeping the assemblies, make keepcontigsonly allows you to do that
keepcontigsonly:
	-find $(ASM_OUT)/* -type f | grep -v $(foreach contigs,$(wildcard $(ALLASMCONTIGS)),-e $(contigs)) -e $(VELVET_OUT_NOSCAF)/noscaf_$(KMAX)/Sequences -e slurm | xargs rm
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
.PRECIOUS: $(VELVETH_OUT_RD) $(VELVETH_OUT_SEQ) $(PRC_READS_OUT)/$(FASTQBASE).qtrim
