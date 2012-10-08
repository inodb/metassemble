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
source $SCRIPTDIR/../../../global-functions.incl

Q1=$1
Q2=$2
QNAME=$3
REF=$4
RNAME=$5
OUTDIR=$6

#TODO: this is some weird uppmax thing that I don't quite follow yet
module load biopython

if [ ! -e $OUTDIR/${RNAME}_${QNAME}-smds.coverage ]; then
    bash $SCRIPTDIR/../../map-bwa-markduplicates.sh -c $Q1 $Q2 $QNAME $REF $RNAME $OUTDIR > /dev/null
fi

python $SCRIPTDIR/gc-content.py $REF | \
awk '{ 
    if (FNR == NR && $1 != "genome") {
        if ($2 == 0) {
            gratio_covered[$1] = 1 - $3 / $4
        }
        gsize[$1] = $4
        gcov[$1] += $5 * $2
    }
    if (FNR != NR) {
        gc[$1] = $2
    } 
}
END {
    print "cov\tlength\tGC_content\tratio_covered"
    for (g in gcov) {
        print g"\t"gcov[g]"\t"gsize[g]"\t"gc[g]"\t"gratio_covered[g]
    }
}' $OUTDIR/${RNAME}_${QNAME}-smds.coverage -
