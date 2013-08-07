#!/bin/bash
HELPDOC=$( cat <<EOF
Scaffold contigs using bambus2. Output-prefix should not contain a path.

Usage:
    bash `basename $0` [options] <input.bam> <contigs.fa> <output-prefix>
Options:
    -k      Keep intermediate files.
    -h      This help documentation.
EOF
)

set -o errexit
set -o nounset
# From: http://tinyurl.com/85qrydz
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPTDIR/../global-functions.incl
COLJAR=$SCRIPTDIR/../../bin/picard-tools-1.77/CollectInsertSizeMetrics.jar
CLEANJAR=$SCRIPTDIR/../../bin/picard-tools-1.77/CleanSam.jar

# Default parameters
RMTMPFILES=true

# Parse options
while getopts ":hk" opt; do
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
shift $(($OPTIND - 1)) 

# Parse arguments
if [ ! $# -eq 3 ]; then
    echo "Invalid number of arguments: 3 needed but $# supplied" >&2
    echo "$HELPDOC"
    exit 1
fi
BAMINPUT=$(readlink -f $1)
BAMBASE=$(basename $BAMINPUT .bam)
BAMDIR=$(dirname $BAMINPUT)
CONTIGS=$(readlink -f $2)
OUTPUTPREFIX=$3
CURDIR=`pwd`
check_prog samtoafg goBambus2 bank-transact samtools
if [ ! -e $COLJAR ]; then
    echo "$COLJAR doesn't exist" >&2
    exit 1
fi
if [ ! -e $CLEANJAR ]; then
    echo "$CLEANJAR doesn't exist" >&2
    exit 1
fi

# TODO: Is it weird that it is in the same dir as the bam, why not an output dir?
cd $BAMDIR

# Convert to sam and filter unmapped reads
samtools view -hF4q1 -o $BAMBASE.sam $BAMBASE.bam
java -jar $CLEANJAR \
    VALIDATION_STRINGENCY=STRICT \
    INPUT=$BAMBASE.sam \
    OUTPUT=$BAMBASE-c.sam 
# Determine insert sizes
java -jar $COLJAR \
    INPUT=$BAMBASE-c.sam \
    OUTPUT=$BAMBASE.insertmetrics \
    HISTOGRAM_FILE=$BAMBASE.inserthist \
    AS=TRUE \
    VALIDATION_STRINGENCY=STRICT 

$SCRIPTDIR/../process-reads/flatten-fasta.sh $CONTIGS > $(basename $CONTIGS).flat
samtoafg \
    -m $(awk '{ if ($5 == "MEAN_INSERT_SIZE") { getline; printf "%i", strtonum($5) } }' $BAMBASE.insertmetrics) \
    -s $(awk '{ if ($6 == "STANDARD_DEVIATION") { getline; printf "%i", strtonum($6) } }' $BAMBASE.insertmetrics) \
    -i 1 -e 1 \
    $(basename $CONTIGS).flat $BAMBASE-c.sam > $BAMBASE.afg
rm -rf $BAMBASE.bnk
bank-transact -cb $BAMBASE.bnk -m $BAMBASE.afg
goBambus2 $BAMBASE.bnk $OUTPUTPREFIX

if $RMTMPFILES; then
    echo "Removing intermediate files."
    rm -rf $BAMBASE.sam \
           $BAMBASE.afg \
           $BAMBASE.bnk \
           $BAMBASE-c.sam \
           $(basename $CONTIGS).flat
fi
