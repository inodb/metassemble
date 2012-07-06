#!/bin/bash
HELPDOC=$( cat <<EOF
Pipeline to run several metagenomic assembly approaches. Currently only accepts
Illumina fastq v1.8 pairs.

= Process reads =
1. Windowed-quality trimming
= De Bruijn Graph Assembly =
2. Multiple velvet assemblies over different kmers with and without -exp_cov auto
= Merging contigs =
3. Merges the assemblies without set exp_cov using cd-hist-est and minimus2

Usage:
    bash `basename $0` [options] <pair1.fastq> <pair2.fastq> <pairbasename>
Options:
    -d      Delete all old output.
    -h      This help documentation.
Example:
    bash `basename $0` fascinating-reads1.fastq fascinating-reads2.fastq fascinating-reads output-folder
EOF
)

set -o errexit
set -o nounset
# From: http://tinyurl.com/85qrydz
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPTDIR/global-functions.incl

function ma_print() {
    echo 'MA:' "$@" | tee -a $LOG
}

# Parse options
DELETEOLD=0
while getopts ":hd" opt; do
    case $opt in
        d)
            DELETEOLD=1
            ;;
        h)
            echo "$HELPDOC"
            exit 0
            ;;
        \?)
            echo "$HELPDOC"
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
    esac
done
shift $(($OPTIND - 1)) 

# Parse arguments
if [ "$#" -ne "4" ]
then
    echo "$HELPDOC"
    echo "Invalid number of arguments: 4 needed but $# supplied" >&2
    exit 1
fi
if [ -f $1 ]; then
    FASTQ1=$1
else
    echo "$HELPDOC"
    echo "$1 is not a file" >&2
    exit 1
fi
if [ -f $2 ]; then
    FASTQ2=$2
else
    echo "$HELPDOC"
    echo "$2 is not a file" >&2
    exit 1
fi
FASTQBASE=$3
OUT=$4
if [ $DELETEOLD -eq 1 ]; then
    rm -rf "$OUT"
fi
PRC_READS_OUT=$OUT/processed-reads
ASM_OUT=$OUT/assemblies
LOG=$OUT/out.log

# Create dirs if they don't exist and check if programs exist
create_dirs $OUT $PRC_READS_OUT $ASM_OUT
check_prog velveth velvetg Trim.pl shuffleSequences_fastq.pl cd-hit-est

# Show input paramaters
ma_print 'Pipeline All'
ma_print 'Pair1 ${FASTQ1}'
ma_print 'Pair2 ${FASTQ2}'
ma_print 'Base name ${FASTQBASE}'
ma_print 'Outfolder ${OUT}'

#echo "Remove duplicates with `which cd-hit-dup`"
#echo "cd-hit-dup -m false -e 0.97 -i ${FASTQ1} -i2 ${FASTQ2}"
#cd-hit-dup -m false -e 0.97 -i ${FASTQ1} -i2 ${FASTQ2}

#  Quality trim
#  qual-type 0 for Sanger quality (or illumina 1.8)
#  type 2 for windowed trimming
if [ ! -f "${PRC_READS_OUT}/${FASTQBASE}.qtrim.unpaired" ] ||
   [ ! -f "${PRC_READS_OUT}/${FASTQ1}.qtrim" ] ||
   [ ! -f "${PRC_READS_OUT}/${FASTQ2}.qtrim" ]; then
ma_print 'Quality trimming'
Trim.pl --qual-type 0 \
        --type 2 \
        --pair1 ${FASTQ1} \
        --pair2 ${FASTQ2} \
        --outpair1 ${PRC_READS_OUT}/${FASTQ1}.qtrim \
        --outpair2 ${PRC_READS_OUT}/${FASTQ2}.qtrim \
        --single ${PRC_READS_OUT}/${FASTQBASE}.qtrim.unpaired \
        | tee -a $LOG
else
    ma_print 'Using previous Trim.pl results'
fi

#  Shuffle sequences to input them into velvet
if [ ! -f ${PRC_READS_OUT}/${FASTQBASE}.qtrim ]; then
    ma_print 'Shuffle sequences for velvet'
    shuffleSequences_fastq.pl ${PRC_READS_OUT}/${FASTQ1}.qtrim ${PRC_READS_OUT}/${FASTQ2}.qtrim ${PRC_READS_OUT}/${FASTQBASE}.qtrim
else
    ma_print 'Using previous shuffleSequences_fastq.pl'
fi

ma_print 'De Bruijn Graph Assembly'
bash $SCRIPTDIR/ma-asm-velvet.sh ${PRC_READS_OUT}/${FASTQBASE}.qtrim ${ASM_OUT}/velvet
#bash $SCRIPTDIR/ma-asm-idba.sh ${PRC_READS_OUT}/${FASTQBASE}.qtrim ${ASM_OUT}/metaidba


# Merge the contigs with minimus2.
#ma_print 'Merge contigs with minimus2'
#create_dirs ${ASM_OUT}/minimus2 ${ASM_OUT}/minimus2/fastq-shortpaired-noscaf-31-40
#bash $SCRIPTDIR/ma-merge-minimus2.sh ${ASM_OUT}/minimus2/fastq-shortpaired-noscaf-31-40 ${ASM_OUT}/velvet/fastq-shortpaired-noscaf_*/contigs.fa
ma_print 'Finished'
