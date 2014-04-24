#!/bin/bash
HELPDOC=$( cat <<EOF
Maps given paired library to given reference with bwa and uses picard to remove
duplicates. BEDTools is used to determine the coverage of the reference.
Returns coverage of reference, GC% content, length and ratio of reference
covered.

Usage:
    bash `basename $0` [options] <reads1> <reads2> <qname> <ref> <rname> <outdir>
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

if [[ ! -e $OUTDIR/${RNAME}_${QNAME}-smds.coverage || $OUTDIR/${RNAME}_${QNAME}-smds.coverage -ot $REF ]]; then
    bash $SCRIPTDIR/../../../map/map-bowtie2-markduplicates.sh -ct 16 $Q1 $Q2 $QNAME $REF $RNAME $OUTDIR > /dev/null
fi

python $SCRIPTDIR/gc-content.py $REF | \
awk \
'{ 
    # go over bedtools histogram file
    if (FNR == NR && $1 != "genome") {
        if ($2 == 0) {
            gratio_covered[$1] = 1 - $3 / $4
        }
        gcov[$1] += $5 * $2
    }
    # go over output from gc-content.py
    if (FNR != NR) {
        nr += 1
        g[nr] = $1
        gc[nr] = $2
        gsize[nr] = $3
    } 
}
END {
    OFS="\t"
    print "cov", "GC_content", "length", "ratio_covered"
    for (i=1; i<= nr; i++) {
        # if there were no reads mapping to the genome, the statistics are not
        # in the bedtools histogram file
        if (!(g[i] in gcov)) {
            gcov[g[i]] = 0
            gratio_covered[g[i]] = 0
        }
        # if there were no bases with 0 coverage in bedtools histogram file,
        # genome was covered completely
        if (!(g[i] in gratio_covered)) {
            gratio_covered[g[i]] = 1
        }
        print g[i], gcov[g[i]], gc[i], gsize[i], gratio_covered[g[i]]
    }
}' $OUTDIR/${RNAME}_${QNAME}-smds.coverage -
