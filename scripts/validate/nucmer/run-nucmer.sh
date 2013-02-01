#!/bin/bash
HELPDOC=$( cat <<EOF
Map two files of contigs against each other using nucmer. The coords file shows
all the alignments, see http://mummer.sourceforge.net/manual/#coords for column
names. The gcoords is the output of show-coords after running delta-filter -g,
which gives you the longest mutually consistent set of alignments for each
ref-query pair. See: http://mummer.sourceforge.net/manual/#filter. 

Usage:
    bash `basename $0` [options] <ref.fa> <query.fa> <output-prefix>
Options:
    -r      Generate report with dnadiff. Report is in <output-prefix>_report.report
    -h      This help documentation.
EOF
)
set -o errexit
set -o nounset
# From: http://tinyurl.com/85qrydz
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPTDIR/../../global-functions.incl

# Default parameters
GEN_REPORT=false

# Parse options
while getopts ":hr" opt; do
    case $opt in
        r)
            GEN_REPORT=true 
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
OUTPUTBASE=$3

check_prog nucmer show-coords

nucmer --maxmatch --prefix=$OUTPUTBASE $REF $QUERY
# delta-filter 
show-coords -rclTH $OUTPUTBASE.delta > $OUTPUTBASE.coords
delta-filter -g $OUTPUTBASE.delta > $OUTPUTBASE.gdelta
show-coords -rclTH $OUTPUTBASE.gdelta > $OUTPUTBASE.gcoords
if $GEN_REPORT; then
    dnadiff -p ${OUTPUTBASE}_report -d $OUTPUTBASE.delta
fi
