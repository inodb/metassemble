#!/bin/bash
HELPDOC=$( cat <<EOF
Maps given paired library to given reference with bwa and uses picard to remove
duplicates. BEDTools is used to determine the coverage of the reference.

Usage:
    bash `basename $0` [options] <reads1> <reads2> <qname> <ref> <rname> <outdir>
Options:
    -k      Keep all output from intermediate steps.
    -h      This help documentation.
EOF
) 

set -o errexit
set -o nounset
# From: http://tinyurl.com/85qrydz
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPTDIR/../global-functions.incl

# Default parameters
RMTMPFILES=true

# Parse options
while getopts "kh" opt; do
    case $opt in
        k)
            RMTMPFILES=false
            ;;
        h)
            echo "$HELPDOC"
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            echo "$HELPDOC"
            exit 1
            ;;
    esac
done

Q1=$1
Q2=$2
QNAME=$3
REF=$4
RNAME=$5
OUTDIR=${6%/}
CURDIR=`pwd`

check_prog bwa samtools genomeCoverageBed

MRKDUP=/bubo/home/h16/inod/glob/downloaded_software/picard-tools-1.66/MarkDuplicates.jar
if [ ! -e $MRKDUP ]; then
    echo "$MRKDUP doesn't exist" >&2
    exit 1
fi

cd $OUTDIR

# Index reference, Burrows-Wheeler Transform
if [ ! -e ${REF}.bwt ]
then
    bwa index $REF
fi

# Find SA coordinates in the reads
bwa aln $REF -1 $Q1 > ${QNAME}1.sai
bwa aln $REF -2 $Q2 > ${QNAME}2.sai

# Align Paired end and bam it
bwa sampe $REF ${QNAME}1.sai ${QNAME}2.sai $Q1 $Q2 > ${RNAME}_${QNAME}.sam
samtools faidx $REF
samtools view -bt $REF.fai ${RNAME}_${QNAME}.sam > ${RNAME}_${QNAME}.bam
samtools sort ${RNAME}_${QNAME}.bam ${RNAME}_${QNAME}-s
samtools index ${RNAME}_${QNAME}-s.bam

# Mark duplicates and sort
java -jar /bubo/home/h16/inod/glob/downloaded_software/picard-tools-1.66/MarkDuplicates.jar \
    INPUT=${RNAME}_${QNAME}-s.bam \
    OUTPUT=${RNAME}_${QNAME}-smd.bam \
    METRICS_FILE=${RNAME}_${QNAME}-smd.metrics \
    AS=TRUE \
    VALIDATION_STRINGENCY=LENIENT \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    REMOVE_DUPLICATES=TRUE
samtools sort ${RNAME}_${QNAME}-smd.bam ${RNAME}_${QNAME}-smds
samtools index ${RNAME}_${QNAME}-smds.bam

# Determine Genome Coverage
genomeCoverageBed -ibam ${RNAME}_${QNAME}-smds.bam > ${RNAME}_${QNAME}-smds.coverage

# Remove temp files
if $RMTMPFILES
then
    rm ${RNAME}_${QNAME}.sam \
       ${RNAME}_${QNAME}.bam \
       ${RNAME}_${QNAME}-smd.bam \
       ${RNAME}_${QNAME}-s.bam \
       ${RNAME}_${QNAME}-s.bam.bai \
       ${QNAME}1.sai \
       ${QNAME}2.sai
fi

cd $CURDIR
