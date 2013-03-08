# Input file arguments.
# This is an example file, change for your library
#
# COMMENTS WITH #
# CAUTION: DON'T WRITE COMMENTS AFTER = SIGN, THEY ARE ASSIGNED TO THE VARIABLE

# Name of input pairs
FASTQ1=pair1.fastq
FASTQ2=pair2.fastq
FASTQBASE=pair

# Output folder
OUT=metassemble

# Set kmin, kmax and stepsize, used for all assemblies: ray velvet
# metavelvet, etc.
KMIN=19
KMAX=75
STEPSIZE=2

# In case input are quality trimmed pairs
# DO_QTRIM=no
# If quality trimmed and interleaved pair is available
# FASTQ_TRIM_IL=pair.fastq

# Create your own target rule
# look in parameters.mk for all possible rules
#
# You can choose for instance from:
# velvet
# metavelvet
# ray
# newbler
# minimus2
# bambus2
subset: velvet metavelvet ray newbler minimus2

# Include metassemble
include $(METASSEMBLE_DIR)/scripts/metassemble.mk
