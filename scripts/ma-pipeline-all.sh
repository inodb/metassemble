#!/bin/bash
HELPDOC=$( cat <<EOF
Pipeline to run several metagenomic assembly approaches. Currently only accepts
Illumina fastq v1.8 pairs.

= Process reads =
Windowed-quality trimming
= De Bruijn Graph Assembly =
Multiple velvet assemblies over different kmers with and without scaffolding.
Multiple metavelvet assemblies over different kmers with and without scaffolding.
Multiple idba assemblies over different kmers with and without scaffolding.
= Merging contigs =
Merge unscaffolded velvet and metavelvet assemblies using cd-hist-est and minimus2
= Validate =
Map all generated assemblies using Mummer and generate plots and statistics on each
assembly.

Usage:
    bash `basename $0` [options] <pair1.fastq> <pair2.fastq> <pairbasename> <output-folder>
Options:
    -d      Delete all old output.
    -h      This help documentation.
Example:
    bash `basename $0` fascinating-reads1.fastq fascinating-reads2.fastq fascinating-reads fascinating-output
EOF
)

set -o errexit
set -o nounset
# From: http://tinyurl.com/85qrydz
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPTDIR/global-functions.incl

function ma_print() {
    echo 'MA:' "$@" | tee -a $LOG
}

# Check if first argument is equal to one of the other
# Return 0 if true, 1 if false
function is_in_array() {
    local step

    for step in "${@:2}"; do
        if [ $1 -eq $step ]; then
            return 0
        fi
    done

    return 1
}

# Parsers a step string of the form 1,2,5-8 and changes it to an array of form 1 2 5 6 7 8
function create_steps_array() {
    if ! [ "$#" -gt 0 ]; then
        echo "No arguments supplied"
    fi

    # Split steps on commas
    steps_spaces=( $(echo $1 | sed 's/,/\ /g') )

    if [[ ${#steps_spaces[@]} -eq 0 ]]; then
        echo "No steps supplied" >&2
        exit 1
    fi

    for step in ${steps_spaces[@]}; do
        if [[ "$step" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
            if ! [ ${steps_all:-unset} = "unset" ]; then
                steps_all=( ${steps_all[@]} $step )
            else
                steps_all=( $step )
            fi
        elif [[ "$step" =~ ^[0-9]+[-][0-9]+$ ]]; then
            low=$(echo $step | cut -d'-' -f1)
            high=$(echo $step | cut -d'-' -f2)
            if ! [[ "$low" =~ ^[0-9]+$ ]] || \
               ! [[ "$high" =~ ^[0-9]+$ ]] || \
                 [ $low -ge $high ]; then
                echo "$step is not of format x-y, where x is a number < y" >&2
                exit 1
            else
                if ! [ ${steps_all:-unset} = "unset" ]; then
                    steps_all=( ${steps_all[@]} $(seq $low $high) )
                else
                    steps_all=( $(seq $low $high) )
                fi
            fi
        else
            echo "$step is not a number or of the form x-y, where x < y" >&2
            exit 1
        fi
    done

    echo ${steps_all[@]}
}

# Parse options
DELETEOLD=false
ALLSTEPS=true
while getopts ":hds:" opt; do
    case $opt in
        d)
            DELETEOLD=true
            ;;
        h)
            echo "$HELPDOC"
            exit 0
            ;;
        s)
            ALLSTEPS=false
            STEPS=( $(create_steps_array $OPTARG) )
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
OUT=$(pwd -P $4)/$4
if $DELETEOLD; then
    rm -rf "$OUT"
fi
PRC_READS_OUT=$OUT/processed-reads
ASM_OUT=$OUT/assemblies
VAL_OUT=$OUT/validate
LOG=$OUT/out.log
#TODO: maybe better to have the validation separate from the assembly process
REF=/bubo/home/h16/inod/metagenomics/results/chris-mock/Project_ID793_dAmore/Sample_50ng_even/mapping/velvet35test/references.fa

# Create dirs if they don't exist and check if programs exist
mkdir -p $OUT $PRC_READS_OUT $ASM_OUT
check_prog velveth velvetg Trim.pl shuffleSequences_fastq.pl cd-hit-est

# Show input paramaters
ma_print 'Pipeline All'
ma_print "Pair1 $FASTQ1"
ma_print "Pair2 $FASTQ2"
ma_print "Base name $FASTQBASE"
ma_print "Output folder $OUT"
if $ALLSTEPS; then
    ma_print "Running all steps"
else
    ma_print "Executing steps: ${STEPS[@]}"
fi

#echo "Remove duplicates with `which cd-hit-dup`"
#echo "cd-hit-dup -m false -e 0.97 -i ${FASTQ1} -i2 ${FASTQ2}"
#cd-hit-dup -m false -e 0.97 -i ${FASTQ1} -i2 ${FASTQ2}

if $ALLSTEPS || $(is_in_array 0 ${STEPS[@]}); then
    #  Quality trim
    #  qual-type 0 for Sanger quality (or illumina 1.8)
    #  type 2 for windowed trimming
    ma_print 'Step 0 Quality trimming'
    Trim.pl --qual-type 0 \
            --type 2 \
            --pair1 ${FASTQ1} \
            --pair2 ${FASTQ2} \
            --outpair1 ${PRC_READS_OUT}/${FASTQBASE}.1.qtrim \
            --outpair2 ${PRC_READS_OUT}/${FASTQBASE}.2.qtrim \
            --single ${PRC_READS_OUT}/${FASTQBASE}.qtrim.unpaired 
fi

if $ALLSTEPS || $(is_in_array 1 ${STEPS[@]}); then
    #  Shuffle sequences (required for e.g. velvet)
    ma_print 'Step 1 shuffle paired sequences'
    shuffleSequences_fastq.pl ${PRC_READS_OUT}/${FASTQBASE}.1.qtrim ${PRC_READS_OUT}/${FASTQBASE}.2.qtrim ${PRC_READS_OUT}/${FASTQBASE}.qtrim
fi
if $ALLSTEPS || $(is_in_array 2 ${STEPS[@]}); then
    ma_print 'Step 2 Velvet Assembly'
    bash $SCRIPTDIR/assembly/run-velvet.sh ${PRC_READS_OUT}/${FASTQBASE}.qtrim ${ASM_OUT}/velvet velveth 31 84 2
    bash $SCRIPTDIR/assembly/run-velvet.sh ${PRC_READS_OUT}/${FASTQBASE}.qtrim ${ASM_OUT}/velvet velvetg 31 84 2
else
    if $(is_in_array 2.1 ${STEPS[@]}); then
        bash $SCRIPTDIR/assembly/run-velvet.sh ${PRC_READS_OUT}/${FASTQBASE}.qtrim ${ASM_OUT}/velvet velveth 31 84 2
    fi
    if $(is_in_array 2.2 ${STEPS[@]}); then
        bash $SCRIPTDIR/assembly/run-velvet.sh ${PRC_READS_OUT}/${FASTQBASE}.qtrim ${ASM_OUT}/velvet velvetg 31 84 2
    fi
fi
if $ALLSTEPS || $(is_in_array 3 ${STEPS[@]}); then
    ma_print 'Step 3 IDBA Assembly'
    bash $SCRIPTDIR/assembly/run-idba.sh ${PRC_READS_OUT}/${FASTQBASE}.qtrim ${ASM_OUT}/idba
fi
if $ALLSTEPS || $(is_in_array 4 ${STEPS[@]}); then
    ma_print 'Step 4 Metavelvet Assembly'
    # TODO: Should maybe check for all velveth output files instead of just KMIN?
    if [ -f ${ASM_OUT}/velvet/fastq-shortpaired-noscaf_31/Sequences ]; then
        mkdir -p ${ASM_OUT}/metavelvet
        cp -r ${ASM_OUT}/velvet/fastq-shortpaired-noscaf_* ${ASM_OUT}/metavelvet/
        bash $SCRIPTDIR/assembly/run-velvet.sh ${PRC_READS_OUT}/${FASTQBASE}.qtrim ${ASM_OUT}/velvet meta-velvetg 31 84 2
    fi
fi
if $ALLSTEPS || $(is_in_array 5 ${STEPS[@]}); then
    ma_print 'Step 5 Merge velvet contigs with minimus2'
    mkdir -p ${ASM_OUT}/minimus2/velvet/noscaf
    bash $SCRIPTDIR/assembly/merge-asm-minimus2.sh ${ASM_OUT}/minimus2/velvet/noscaf ${ASM_OUT}/velvet/*noscaf*/contigs.fa
fi
if $ALLSTEPS || $(is_in_array 6 ${STEPS[@]}); then
    ma_print 'Step 6 Merge metavelvet contigs with minimus2'
    mkdir -p ${ASM_OUT}/minimus2/metavelvet/noscaf
    bash $SCRIPTDIR/assembly/merge-asm-minimus2.sh ${ASM_OUT}/minimus2/metavelvet/noscaf ${ASM_OUT}/metavelvet/*noscaf*/contigs.fa
fi
if $ALLSTEPS || $(is_in_array 7 ${STEPS[@]}); then
    ma_print 'Step 7 Determine coverage and GC content reference'
    mkdir -p ${VAL_OUT}/reference
    # TODO: Tidy up
    if [ ! -f $VAL_OUT/reference/ref.stats ]; then
        #if $SBATCH; then
        #    prefix="sbatch -A b2010008 -J map.ref.$FASTQBASE -t 1-00:00:00 -p node /bubo/home/h16/inod/bin/sbatch_job"
        #else
        #    prefix='eval'
        #fi
        #$prefix \
        bash $SCRIPTDIR/validate/reference/stats/length-gc-cov.sh ${PRC_READS_OUT}/${FASTQBASE}.1.qtrim ${PRC_READS_OUT}/${FASTQBASE}.2.qtrim ${FASTQBASE}.qtrim $REF ref $VAL_OUT/reference > $VAL_OUT/reference/ref.stats
    else
        ma_print "$VAL_OUT/reference/ref.stats exists"        
    fi
fi

function validate() {
    # Validate given assemblies.
    # First argument is the root output dir, each subfolder gets the same name as the folder that the assembly is in.
    for asm in ${@:2}; do
        local out
        out=$1/$(basename $(dirname $asm))
        mkdir -p $out
        
        #if $SBATCH; then
        #    prefix="sbatch -A b2010008 -J val.$(basename $(dirname $(dirname $asm))).$FASTQBASE -t 01:00:00 -p node /bubo/home/h16/inod/bin/sbatch_job"
        #else
        #    prefix='eval'
        #fi

        # Only run mummer if it has not been run before
        if [ ! -f $out/nucmer.gcoords ]; then
            bash $SCRIPTDIR/validate/nucmer/run-nucmer.sh $REF $asm $out/nucmer;
        fi;
        bash $SCRIPTDIR/validate/nucmer/calc-stats-n-plot.sh $REF $asm $out/nucmer.gcoords $VAL_OUT/reference/ref.stats $out
    done
}

if $ALLSTEPS || $(is_in_array 8 ${STEPS[@]}); then
    ma_print 'Step 8 Validate velvet'
    validate $VAL_OUT/assemblies/velvet ${ASM_OUT}/velvet/*/contigs.fa
fi
if $ALLSTEPS || $(is_in_array 9 ${STEPS[@]}); then
    ma_print 'Step 9 Validate metavelvet'
    validate $VAL_OUT/assemblies/metavelvet ${ASM_OUT}/metavelvet/*/contigs.fa
fi
if $ALLSTEPS || $(is_in_array 10 ${STEPS[@]}); then
    ma_print 'Step 10 Validate idba'
    validate $VAL_OUT/assemblies/idba ${ASM_OUT}/idba/*/*-contig.fa
fi
if $ALLSTEPS || $(is_in_array 11 ${STEPS[@]}); then
    ma_print 'Step 11 Validate minimus2 velvet merged contigs'
    validate $VAL_OUT/assemblies/minimus2/velvet ${ASM_OUT}/minimus2/velvet/*/all-merged.fasta
fi
if $ALLSTEPS || $(is_in_array 12 ${STEPS[@]}); then
    ma_print 'Step 12 Validate minimus2 metavelvet merged contigs'
    validate $VAL_OUT/assemblies/minimus2/metavelvet ${ASM_OUT}/minimus2/metavelvet/*/all-merged.fasta
fi
if $ALLSTEPS || $(is_in_array 13 ${STEPS[@]}); then
    ma_print 'Step 13 Plot purity versus contig length for all assemblies'
    #TODO Everything is using names hard-coded in run-velvet.sh, this should change.
    (
        awk '{ if (NR == 2) { print $0,"kmer_size","kmer_type" } }' $VAL_OUT/assemblies/velvet/fastq-shortpaired-noscaf_31/purity;
        for i in {31..83..2}; do 
            if [ -f "$VAL_OUT/assemblies/metavelvet/fastq-shortpaired-noscaf_${i}/purity" ]; then
                awk '{ if (NR == 3) { $1="metavelvetnoscaf"; print $0,'$i',"single" } }' \
                    $VAL_OUT/assemblies/metavelvet/fastq-shortpaired-noscaf_${i}/purity;
            fi
        done;
        for i in {31..83..2}; do 
            if [ -f "$VAL_OUT/assemblies/metavelvet/fastq-shortpaired-scaf_${i}/purity" ]; then
                awk '{ if (NR == 3) { $1="metavelvetscaf"; print $0,'$i',"single" } }' \
                    $VAL_OUT/assemblies/metavelvet/fastq-shortpaired-scaf_${i}/purity;
            fi
        done;
        for i in {31..83..2}; do 
            if [ -f "$VAL_OUT/assemblies/velvet/fastq-shortpaired-noscaf_${i}/purity" ]; then
                awk '{ if (NR == 3) { $1="velvetnoscaf"; print $0,'$i',"single" } }' \
                    $VAL_OUT/assemblies/velvet/fastq-shortpaired-noscaf_${i}/purity;
            fi
        done;
        for i in {31..83..2}; do
            
            if [ -f "$VAL_OUT/assemblies/velvet/fastq-shortpaired-scaf_${i}/purity" ]; then
                awk '{ if (NR == 3) { $1="velvetscaf"; print $0,'$i',"single" } }' \
                    $VAL_OUT/assemblies/velvet/fastq-shortpaired-scaf_${i}/purity;
            fi
        done;
        awk '{ if (NR == 3) { $1="minimus2velvetnoscaf"; print $0,"31","multiple";  print $0,"83","multiple" } }' \
            $VAL_OUT/assemblies/minimus2/velvet/noscaf/purity;
        awk '{ if (NR == 3) { $1="minimus2metavelvetnoscaf"; print $0,"31","multiple";  print $0,"83","multiple" } }' \
            $VAL_OUT/assemblies/minimus2/metavelvet/noscaf/purity;
    ) >  $VAL_OUT/assemblies/purity
    Rscript $SCRIPTDIR/validate/nucmer/plot/purity-L50.R --gradient-color-column kmer_size --title $FASTQBASE $VAL_OUT/assemblies/l50-purity-$FASTQBASE.pdf $VAL_OUT/assemblies/purity
fi

ma_print 'Finished'
