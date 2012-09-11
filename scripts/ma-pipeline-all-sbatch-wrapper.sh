#!/bin/bash
HELPDOC=$( cat <<EOF
Schedules assembly steps from ma-pipelin-all.sh using sbatch for parallel
computation. Change script according to your needs for the sequences your
working with.

Usage:
    bash `basename $0` [options] <pair1.fastq> <pair2.fastq> <pairbasename> <output-folder>
EOF
)

set -o errexit
set -o nounset
# From: http://tinyurl.com/85qrydz
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Parse options
while getopts ":h" opt; do
    case $opt in
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
if [ "$#" -ne "4" ]; then
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
OUT=$(cd $4 && pwd)

function ma_pipeline_wrapper() {
    sbatch -J $1.$FASTQBASE --output=slurm-$1.out --error=slurm-$1.err -A b2010008 $2 /bubo/home/h16/inod/bin/sbatch_job bash $SCRIPTDIR/ma-pipeline-all.sh -s $3 $FASTQ1 $FASTQ2 $FASTQBASE $OUT | cut -d' ' -f4
}

set -o xtrace
trimmingjobid="$(ma_pipeline_wrapper qtrim "-t 4:00:00 -p core" 0,1)"
velvethjobid=$(ma_pipeline_wrapper velveth "-t 1-00:00:00 -p node -C mem24GB -d afterok:$trimmingjobid" 2.1)
velvetgjobid=$(ma_pipeline_wrapper velvetg "-t 03:00:00 -p node -d afterok:$velvethjobid" 2.2)
metavelvetgjobid=$(ma_pipeline_wrapper metavelvetg "-t 03:00:00 -p node -d afterok:$velvethjobid" 4)
minimus2velvetjobid=$(ma_pipeline_wrapper minimus2velvet "-t 04:00:00 -p node -d afterok:$velvetgjobid" 5)
minimus2metavelvetjobid=$(ma_pipeline_wrapper minimus2metavelvet "-t 04:00:00 -p node -d afterok:$velvetgjobid" 6)
set +o xtrace
