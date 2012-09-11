set -o errexit
set -o nounset
# From: http://tinyurl.com/85qrydz
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPTDIR/../global-functions.incl

if [ -f $1 ]; then
    FASTQPAIRED=$1
else
    #echo "$HELPDOC"
    echo "$1 is not a file" >&2
    exit 1
fi
OUTDIR=${2%/}

check_prog idba

function run_idba() {
    local out
    parameters=`echo "${@}" | tr -d ' '`
    create_dirs $OUTDIR/$parameters
    
    set -o xtrace
    idba \
      --read $FASTQPAIRED \
      --output $OUTDIR/$parameters/out \
      $@
    set +o xtrace
}

create_dirs $OUTDIR
run_idba --mink 25 --maxk 45
run_idba --mink 45 --maxk 65
run_idba --mink 65 --maxk 75
run_idba --mink 25 --maxk 75
run_idba --mink 25 --maxk 45 --scaffold
run_idba --mink 45 --maxk 65 --scaffold
run_idba --mink 65 --maxk 75 --scaffold
run_idba --mink 25 --maxk 75 --scaffold

# Remove binary kmer files, because they are huge and unreadable anyway
rm $OUTDIR/*/*.kmer
