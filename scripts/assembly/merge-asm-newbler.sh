#!/bin/bash
HELPDOC=$( cat <<EOF
Merges contigs using Newbler.

Usage:
    bash `basename $0` [options] <output directory> <file1.fasta> <file2.fasta> [filen.fasta ...]
Options:
    -h      This help documentation.
EOF
)
set -o errexit
set -o nounset
# From: http://tinyurl.com/85qrydz
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
source $SCRIPTDIR/global-functions.incl

# Parse options
while getopts ":h" opt; do
    case $opt in
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
if [ "$#" -lt "3" ]
then
    echo "Invalid number of arguments: 3 needed but $# supplied" >&2
    echo "$HELPDOC"
    exit 1
fi
if [ ! -d $1 ]; then
    mkdir $1
fi
OUTDIR=${1%/}
CURDIR=`pwd`

create_dirs $OUTDIR

module load biopython
python $SCRIPTDIR/process-reads/cut-up-fasta.py ${@:2} > $OUTDIR/cut-up.fasta
module load 454-dataanalysis/2.6
runAssembly -force -o $OUTDIR $OUTDIR/cut-up.fasta
rm $OUTDIR/cut-up.fasta
