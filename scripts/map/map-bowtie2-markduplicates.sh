#!/bin/bash
HELPDOC=$( cat <<EOF
Maps given paired library to given reference with bowtie2 and uses picard to remove
duplicates.

Usage:
    bash `basename $0` [options] <reads1> <reads2> <qname> <ref> <rname> <outdir>
        
        <reads1> fast{q,a,q.gz,a.gz} file with R1 reads
        <reads2> fast{q,a,q.gz,a.gz} file with R2 reads
        <qname>  name of reads - used for output name
        <ref>    fasta file with reference fasta
        <rname>  name of ref - used for output name
Options:
    -t      Number of threads for bowtie2 and the java garbage collector
    -c      Calculate coverage with BEDTools
    -k      Keep all output from intermediate steps.
    -h      This help documentation.
    -p      Set mapping parameters for bowtie2
            (default: nothing i.e. bowtie default parameters)
    -j      Set java parameters for MarkDuplicates
            (default: -Xms2g -Xmx32g -XX:ParallelGCThreads=NR_THREADS -XX:MaxPermSize=2g -XX:+CMSClassUnloadingEnabled)
Example:
    map-bowtie2-markduplicates.sh -ct 16 reads_R1.fastq.gz reads_R2.fastq.gz pair contigs.fa asm bowtie2
Output:
    All output is in <outdir>:

    <rname>_<qname>-smd.metrics             - the metrics output from removing
                                              PCR duplicates with
                                              MarkDuplicates (smd stands for,
                                              samtools sort, followed by
                                              MarkDuplicates)
    <rname>_<qname>-smds.bam                - the bam output from mapping the
                                              reads and removing the duplicates
                                              with MarkDuplicates, followed by
                                              one more time sorting
    <rname>_<qname>-smds.flagstat           - results form samtools flagstat on
                                              bam file, gives read mapping stats

    If -c is specified (requires genomeCoverageBed and BioPython):

    <rname>_<qname>-smds.coverage           - result from running
                                              genomeCoverageBed -ibam on the
                                              bam file, gives histogram
                                              coverage per contig
    <rname>_<qname>-smds.coverage.percontig - gives mean coverage per contig
    
EOF
) 

set -o errexit
set -o nounset
# From: http://tinyurl.com/85qrydz
#SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#source $SCRIPTDIR/../global-functions.incl
#MRKDUP=$SCRIPTDIR/../../bin/picard-tools-1.77/MarkDuplicates.jar
source $METASSEMBLE_DIR/scripts/global-functions.incl
MRKDUP=$METASSEMBLE_DIR/bin/picard-tools-1.77/MarkDuplicates.jar

# Default parameters
RMTMPFILES=true
CALCCOV=false
THREADS=1
BOWTIE2_OPT=''
JAVA_OPT=

# Parse options
while getopts "khct:p:j:" opt; do
    case $opt in
        c)
            CALCCOV=true
            ;;
        k)
            RMTMPFILES=false
            ;;
        t)
            THREADS=$OPTARG
            ;;
        p)
            BOWTIE2_OPT=$OPTARG
            ;;
        j)
            JAVA_OPT=$OPTARG
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

if [[ -z $JAVA_OPT ]]
then
    JAVA_OPT="-Xms2g -Xmx32g -XX:ParallelGCThreads=$THREADS -XX:MaxPermSize=2g -XX:+CMSClassUnloadingEnabled"
fi

if [ "$#" -ne "6" ]
then
    echo "Invalid number of arguments: 6 needed but $# supplied" >&2
    echo "$HELPDOC"
    exit 1
fi
Q1=$1
Q2=$2
QNAME=$3
REF=$4
RNAME=$5
OUTDIR=${6%/}
#CURDIR=`pwd`

check_prog bowtie2 samtools genomeCoverageBed

if [ ! -e $MRKDUP ]; then
    echo "$MRKDUP doesn't exist" >&2
    exit 1
fi

mkdir -p $OUTDIR

# Index reference, Burrows-Wheeler Transform
if [ ! -e ${REF}.1.bt2 ]
then
    bowtie2-build $REF $REF
fi

# Align Paired end and bam it
bowtie2 ${BOWTIE2_OPT} -p $THREADS -x $REF -1 $Q1 -2 $Q2 -S $OUTDIR/${RNAME}_${QNAME}.sam
# Index reference for samtools
if [ ! -e ${REF}.fai ]
then
    samtools faidx $REF
fi
samtools view -bt $REF.fai $OUTDIR/${RNAME}_${QNAME}.sam > $OUTDIR/${RNAME}_${QNAME}.bam
samtools sort $OUTDIR/${RNAME}_${QNAME}.bam $OUTDIR/${RNAME}_${QNAME}-s
samtools index $OUTDIR/${RNAME}_${QNAME}-s.bam

# Mark duplicates and sort
java $JAVA_OPT \
    -jar $MRKDUP \
    INPUT=$OUTDIR/${RNAME}_${QNAME}-s.bam \
    OUTPUT=$OUTDIR/${RNAME}_${QNAME}-smd.bam \
    METRICS_FILE=$OUTDIR/${RNAME}_${QNAME}-smd.metrics \
    AS=TRUE \
    VALIDATION_STRINGENCY=LENIENT \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    REMOVE_DUPLICATES=TRUE
samtools sort $OUTDIR/${RNAME}_${QNAME}-smd.bam $OUTDIR/${RNAME}_${QNAME}-smds
samtools index $OUTDIR/${RNAME}_${QNAME}-smds.bam
samtools flagstat $OUTDIR/${RNAME}_${QNAME}-smds.bam > $OUTDIR/${RNAME}_${QNAME}-smds.flagstat

# Determine Genome Coverage and mean coverage per contig
if $CALCCOV; then
    genomeCoverageBed -ibam $OUTDIR/${RNAME}_${QNAME}-smds.bam > $OUTDIR/${RNAME}_${QNAME}-smds.coverage
    # generate table with length and coverage stats per contig (From http://github.com/BinPro/CONCOCT)
    python $METASSEMBLE_DIR/scripts/validate/map/gen_contig_cov_per_bam_table.py --isbedfiles $REF $OUTDIR/${RNAME}_${QNAME}-smds.coverage > $OUTDIR/${RNAME}_${QNAME}-smds.coverage.percontig
fi

# Remove temp files
if $RMTMPFILES; then
    rm $OUTDIR/${RNAME}_${QNAME}.sam \
       $OUTDIR/${RNAME}_${QNAME}.bam \
       $OUTDIR/${RNAME}_${QNAME}-smd.bam \
       $OUTDIR/${RNAME}_${QNAME}-s.bam \
       $OUTDIR/${RNAME}_${QNAME}-s.bam.bai
fi
