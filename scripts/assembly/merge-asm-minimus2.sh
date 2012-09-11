#!/bin/bash
HELPDOC=$( cat <<EOF
Merges contigs using minimus2 as in http://ged.msu.edu/angus/metag-assembly-2011/velvet-multik.html.

Usage:
    bash `basename $0` [options] <output directory> <file1.fasta> <file2.fasta> [filen.fasta ...]
Options:
    -h      This help documentation.
EOF
)

set -o errexit
set -o nounset
# From: http://tinyurl.com/85qrydz
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPTDIR/../global-functions.incl

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

check_prog toAmos minimus2 cd-hit-est
create_dirs $OUTDIR

echo 'Concatenate contigs'
cat ${@:2} > $OUTDIR/concatenated.fasta
echo `grep -c '>' $OUTDIR/concatenated.fasta` 'contigs'
echo 'Done'

# Rename the contigs (have to be unique names) and only keep contigs >= 200bp
cd $OUTDIR
echo 'Rename contigs (to make sure they are unique)'
#TODO: this is some weird uppmax thing that I don't quite follow yet
module load biopython
python $SCRIPTDIR/../process-reads/reads-rename.py -oc 200 concatenated.fasta > concatenated-min200-rename.fasta
echo `grep -c '>' concatenated.fasta` 'contigs'
echo 'Done'

# Remove 99% similar sequences else minimus2 has some issues merging contigs
# (Might be something like it doesn't know which one to merge with cause
# multiple contigs are so similar)
echo 'Remove similar sequences'
cd-hit-est -i concatenated-min200-rename.fasta \
           -o concatenated-min200-rename-derep.fasta \
           -c 0.99 \
           -M 7500
echo 'Done'

# Merge contigs with minimus2
echo 'Merge contigs with minimus2'
toAmos -s concatenated-min200-rename-derep.fasta \
       -o concatenated-min200-rename-derep.afg
minimus2 concatenated-min200-rename-derep
# Combine merged and non-merged contigs
cat concatenated-min200-rename-derep.fasta \
    concatenated-min200-rename-derep.singletons.seq > all-merged.fasta
cd $CURDIR
echo 'Done'
