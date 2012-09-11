set -o errexit
set -o nounset
# From: http://tinyurl.com/85qrydz
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPTDIR/../../global-functions.incl

if [[ -f $1 ]]; then
    REF=$1
else
    echo "$1 is not a file" >&2
    exit 1
fi
if [[ -f $2 ]]; then
    QUERY=$2
else
    echo "$2 is not a file" >&2
    exit 1
fi
if [[ -f $3 ]]; then
    COORDS=$3
else
    echo "$3 is not a file" >&2
    exit 1
fi
if [[ -f $4 ]]; then
    REFSTATS=$4
else
    echo "$4 is not a file" >&2
    exit 1
fi
OUTDIR=${5%/}

set -o xtrace 2>/dev/null
bash $SCRIPTDIR/process-reads/inodb-flatten-fasta.sh $QUERY | gawk -f $SCRIPTDIR/stats/purity.awk - $COORDS > $OUTDIR/purity
Rscript $SCRIPTDIR/stats/pergenome-perassembly.R $QUERY $REFSTATS $COORDS $OUTDIR
set +o xtrace 2>/dev/null
