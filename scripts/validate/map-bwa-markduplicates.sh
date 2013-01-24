#!/bin/bash
HELPDOC=$( cat <<EOF
Maps given paired library to given reference with bwa and uses picard to remove
duplicates.

Usage:
    bash `basename $0` [options] <reads1> <reads2> <qname> <ref> <rname> <outdir>
Options:
    -a      Under development: bwa aln options seperated by comma e.g. -a -o0,-e0,-n0,-k0
            for exact matching
    -t      Number of threads for bwa aln
    -c      Calculate coverage with BEDTools
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
CALCCOV=false
THREADS=1
BWA_ALN_OPT=""

# Parse options
while getopts "a:khct:" opt; do
    case $opt in
        a)
            #TODO: Proper bwa aln option catching?
            #Example exact aln: -o0,-e0,-n0,-k0
            BWA_ALN_OPT=`echo $OPTARG | sed 's/,/\ /g'`
            ;;
        c)
            CALCCOV=true
            ;;
        k)
            RMTMPFILES=false
            ;;
        t)
            THREADS=$OPTARG
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
shift $(($OPTIND - 1)) 

if [ "$#" -ne "6" ]
then
    echo "Invalid number of arguments: 6 needed but $# supplied" >&2
    echo "$HELPDOC"
    exit 1
fi
Q1=$(readlink -f $1)
if [ ! -f "$Q1" ]
then
    echo "Pair 1 doesn't exist: $1"
    exit 1
fi
Q2=$(readlink -f $2)
if [ ! -f "$Q2" ]
then
    echo "Pair 2 doesn't exist: $2"
    exit 1
fi
QNAME=$3
REF=$(readlink -f $4)
if [ ! -f "$REF" ]
then
    echo "Reference doesn't exist: $4"
    exit 1
fi
RNAME=$5
OUTDIR=${6%/}
CURDIR=`pwd`

check_prog bwa samtools genomeCoverageBed samtoafg

MRKDUP=/bubo/home/h16/inod/glob/downloaded_software/picard-tools-1.77/MarkDuplicates.jar
if [ ! -e $MRKDUP ]; then
    echo "$MRKDUP doesn't exist" >&2
    exit 1
fi

mkdir -p $OUTDIR
cd $OUTDIR

# Index reference, Burrows-Wheeler Transform
if [ ! -e ${REF}.bwt ]
then
    bwa index $REF
fi

if [[ ! -s $Q1 || ! -s $Q2 ]]; then
    echo "$Q1 or $Q2 is empty" >&2
    exit 1
fi

# Find SA coordinates in the reads
bwa aln $BWA_ALN_OPT -t $THREADS $REF -1 $Q1 > ${QNAME}1.sai
bwa aln $BWA_ALN_OPT -t $THREADS $REF -2 $Q2 > ${QNAME}2.sai

# Align Paired end and bam it
bwa sampe $REF ${QNAME}1.sai ${QNAME}2.sai $Q1 $Q2 > ${RNAME}_${QNAME}.sam
samtools faidx $REF
samtools view -bt $REF.fai ${RNAME}_${QNAME}.sam > ${RNAME}_${QNAME}.bam
samtools sort ${RNAME}_${QNAME}.bam ${RNAME}_${QNAME}-s
samtools index ${RNAME}_${QNAME}-s.bam

# Mark duplicates and sort
java -jar $MRKDUP \
    INPUT=${RNAME}_${QNAME}-s.bam \
    OUTPUT=${RNAME}_${QNAME}-smd.bam \
    METRICS_FILE=${RNAME}_${QNAME}-smd.metrics \
    AS=TRUE \
    VALIDATION_STRINGENCY=LENIENT \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    REMOVE_DUPLICATES=TRUE
samtools sort ${RNAME}_${QNAME}-smd.bam ${RNAME}_${QNAME}-smds
samtools index ${RNAME}_${QNAME}-smds.bam

# Determine Genome Coverage and mean coverage per contig
if $CALCCOV; then
    genomeCoverageBed -ibam ${RNAME}_${QNAME}-smds.bam > ${RNAME}_${QNAME}-smds.coverage
    awk 'BEGIN {pc=""} 
    {
        c=$1;
        if (c == pc) {
            cov=cov+$2*$5;
        } else {
            print pc,cov;
            cov=$2*$5;
        pc=c}
    } END {print pc,cov}' ${RNAME}_${QNAME}-smds.coverage | tail -n +2 > ${RNAME}_${QNAME}-smds.coverage.percontig
fi

# Remove temp files
if $RMTMPFILES; then
    rm ${RNAME}_${QNAME}.sam \
       ${RNAME}_${QNAME}.bam \
       ${RNAME}_${QNAME}-smd.bam \
       ${RNAME}_${QNAME}-s.bam \
       ${RNAME}_${QNAME}-s.bam.bai \
       ${QNAME}1.sai \
       ${QNAME}2.sai
fi

cd $CURDIR
