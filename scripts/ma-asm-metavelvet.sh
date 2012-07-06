set -o errexit
set -o nounset
# From: http://tinyurl.com/85qrydz
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPTDIR/global-functions.incl

if [ -d $1 ]; then
    VELVETASM=$1
else
    #echo "$HELPDOC"
    echo "$1 is not a directory" >&2
    exit 1
fi

check_prog meta-velvetg
meta-velvetg $VELVETASM
