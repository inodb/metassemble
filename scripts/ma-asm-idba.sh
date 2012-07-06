set -o errexit
set -o nounset
# From: http://tinyurl.com/85qrydz
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPTDIR/global-functions.incl

if [ -f $1 ]; then
    FASTQPAIRED=$1
else
    #echo "$HELPDOC"
    echo "$1 is not a file" >&2
    exit 1
fi
OUTDIR=${2%/}
KMIN=31
KMAX=84

MINCOUNTMIN=2
MINCOUNTMAX=20
MINCOUNTSTEP=2

MINPAIRMIN=2
MINPAIRMAX=10
MINPAIRSTEP=2

check_prog idba

create_dirs $OUTDIR
CURDIR=`pwd`
for ((i=$MINPAIRMIN; $i<$MINPAIRMAX ;i=$i+$MINPAIRSTEP))
do
    for ((j=$MINCOUNTMIN; $j<$MINCOUNTMAX; j=$j+$MINCOUNTSTEP))
    do
        create_dirs $OUTDIR/minCount${j}minPairs${i}
        idba \
          --read $FASTQPAIRED \
          --output $OUTDIR/minCount${j}minPairs${i}/minCount${j}minPairs${i} \
          --mink $KMIN \
          --maxk $KMAX \
          --minPairs $i \
          --minCount $j 
    done
done
