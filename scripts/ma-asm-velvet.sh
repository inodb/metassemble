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
STEPSIZE=2
check_prog velveth velvetg

# Velveth over multiple kmers
create_dirs $OUTDIR
velveth $OUTDIR/fastq-shortpaired-noscaf \
       $KMIN,$KMAX,$STEPSIZE \
       -fastq -shortPaired ${FASTQPAIRED}
# Velvet has three primary parameters, the k size, exp_cov and cov_cutoff.
# See also:
# http://www.homolog.us/blogs/2012/06/08/an-explanation-of-velvet-parameter-exp_cov/
# Since velvet wasn't made for metagenomics we will not use these, we do need
# the exp_cov auto results for meta-velvetg
# velvetg no scaffolding, no exp_cov and exp_cov auto
i=$KMIN
for ((i=$KMIN; $i<$KMAX; i=$i+$STEPSIZE))
do
    # Copy velveth results
    cp -r ${OUTDIR}/fastq-shortpaired-noscaf_${i} \
          ${OUTDIR}/fastq-shortpaired-noscaf-expcovauto_${i}
    # no scaf, no exp_cov
    velvetg $OUTDIR/fastq-shortpaired-noscaf_${i} \
        -scaffolding no
    # no scaf, exp_cov auto
    velvetg ${OUTDIR}/fastq-shortpaired-noscaf-expcovauto_${i} \
        -scaffolding no \
        -exp_cov auto
    # Run meta-velvet
    bash $SCRIPTDIR/ma-asm-metavelvet.sh \
        ${OUTDIR}/fastq-shortpaired-noscaf-expcovauto_${i}
done
