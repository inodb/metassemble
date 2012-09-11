set -o errexit
set -o nounset
# From: http://tinyurl.com/85qrydz
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPTDIR/../global-functions.incl
HELPDOC=$( cat <<EOF
Runs velveth, velvetg or meta-velvetg on paired reads with default parameters
over KMIN to KAX with given STEPSZIE with and without scaffolding. Results are
in fastq-shortpaired-noscaf for no scaffolding and in fastq-shortpaired-scaf for
scaffolding. Meta-velvetg and velvetg can easily be run parallel. The different
velvetg and meta-velvetg assemblies can be run parallel as well by calling the
script with equal KMIN and KMAX values \(STEPSIZE doesn\'t matter in that
case\). In the latter case, you might as well just run velvet without this
wrapper script though. The entire existence of this script is arguably
useless, but so is yours.

Usage:
    bash `basename $0` [options] <interleaved-pairs.fastq> <output-folder> <(velveth|velvetg|meta-velvetg)> KMIN KMAX STEPSIZE
EOF
)

# Parse arguments
if [ "$#" -ne "6" ]; then
    echo "$HELPDOC"
    echo "Invalid number of arguments: 6 needed but $# supplied" >&2
    exit 1
fi
if [ -f $1 ]; then
    FASTQPAIRED=$1
else
    echo "$HELPDOC"
    echo "$1 is not a file" >&2
    exit 1
fi
OUTDIR=${2%/}
PROGRAM=$3
if [ ! "$PROGRAM" = "velveth" ] && [ ! "$PROGRAM" = "velvetg" ] && [ ! "$PROGRAM" = "meta-velvetg" ]; then
    echo "$PROGRAM is invalid program specification, must be velvetg or velveth" >&2
    exit 1
fi
KMIN=$4
KMAX=$5
STEPSIZE=$6

# Check if programs are available
check_prog velveth velvetg meta-velvetg

# Velveth over multiple kmers
create_dirs $OUTDIR
# If sbatch, run as sbatch job
if [ "$PROGRAM" = "velveth" ]; then
    velveth $OUTDIR/fastq-shortpaired-noscaf \
           $KMIN,$KMAX,$STEPSIZE \
           -fastq -shortPaired ${FASTQPAIRED}
fi
# Velvet has three primary parameters, the k size, exp_cov and cov_cutoff.
# See also:
# http://www.homolog.us/blogs/2012/06/08/an-explanation-of-velvet-parameter-exp_cov/
# Since velvet wasn't made for metagenomics we will not use any.
# Run over several kmers with and without scaffolding
if $SBATCH; then
    prefix="sbatch --output=slurm-velvetg.output --error=slurm-velvetg.err -A b2010008 -J velvetg.$(basename $FASTQPAIRED) -t 02:00:00 -p node -d afterok:$velvethjobid /bubo/home/h16/inod/bin/sbatch_job "
else
    prefix='eval'
fi
# no scaf and scaf for velvetg
if [ "$PROGRAM" = "velvetg" ]; then
    for asm in $OUTDIR/fastq-shortpaired-noscaf_*; do
        velvetg \$asm -scaffolding no || { echo \$asm -scaffolding no failed. >&2; rm -rf \$asm; continue; };
        scafasm=\$(echo \$asm | sed 's/noscaf/scaf/');
        cp -r \$asm \$scafasm;
        velvetg \$scafasm -scaffolding yes || { echo \$scafasm -scaffolding yes failed. >&2; rm -rf \$scafasm; continue; };
    done
fi
# no scaf and scaf for meta-velvetg (note this overwrites possible previous velvetg contigs)
if [ "$PROGRAM" = "meta-velvetg" ]; then
    for asm in $OUTDIR/fastq-shortpaired-noscaf_*; do
        velvetg \$asm -scaffolding no -exp_cov auto -read_trkg yes || { echo \$asm -scaffolding no -exp_cov auto -read_trkg yes failed. >&2; rm -rf \$asm; exit 1; };
        meta-velvetg \$asm -scaffolding no || { echo \$asm -scaffolding no failed. >&2; rm -rf \$asm; continue; };
        scafasm=\$(echo \$asm | sed 's/noscaf/scaf/');
        cp -r \$asm \$scafasm;
        meta-velvetg \$scafasm -scaffolding yes || { echo \$scafasm -scaffolding yes failed. >&2; rm -rf \$scafasm; continue; };
    done
fi
