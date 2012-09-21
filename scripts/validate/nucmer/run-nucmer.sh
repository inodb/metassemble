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
OUTPUTBASE=$3

check_prog nucmer show-coords

nucmer --maxmatch --prefix=$OUTPUTBASE $REF $QUERY
# delta-filter 
show-coords -rclTH $OUTPUTBASE.delta > $OUTPUTBASE.coords
delta-filter -g $OUTPUTBASE.delta > $OUTPUTBASE.gdelta
show-coords -rclTH $OUTPUTBASE.gdelta > $OUTPUTBASE.gcoords
