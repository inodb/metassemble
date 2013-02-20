#!/bin/bash
HELPDOC=$( cat <<EOF
Checks if all programs necessary to run MetAssemble are installed.

Usage:
    bash `basename $0`
EOF
) 
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
if [ "$#" -gt "0" ]
then
    echo "Invalid number of arguments: 0 needed but $# supplied" >&2
    echo "$HELPDOC"
    exit 1
fi

#set -o errexit
set -o nounset
# From: http://tinyurl.com/85qrydz
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Check if programs are installed
function check_prog_verbose() {
    local progloc
    progloc=`which $1 2>/dev/null`
    if [ $? -ne 0 ]; then
        echo $2 >&2
        exit 1
    else
        echo "Using: $progloc"
    fi
}

# Check if programs are installed
function check_prog_python_verbose() {
    local progloc
    python -c "import $1" 2>/dev/null
    if [ $? -ne 0 ]; then
        echo $2 >&2
        exit 1
    else
        echo "Python module $1 installed"
    fi
}

function test_Ray() {
    local output
    output=`Ray --help 2>&1`
    if [ $? -ne 0 ]; then
        echo "Ray test failed, see output:" >&2
        echo $output >&2
        exit 1
    else
        echo "Ray seems to be working"
    fi
}

echo "Check dependencies..."
check_prog_verbose sickle "sickle missing from PATH, download from https://github.com/najoshi/sickle"
check_prog_verbose velvetg "velvetg missing from PATH, download from http://www.ebi.ac.uk/~zerbino/velvet"
check_prog_verbose velveth "velveth missing from PATH, download from http://www.ebi.ac.uk/~zerbino/velvet"
check_prog_verbose meta-velvetg "meta-velvetg missing from PATH, download from http://metavelvet.dna.bio.keio.ac.jp/"
check_prog_verbose Ray "Ray missing from PATH, download from http://denovoassembler.sourceforge.net/"
test_Ray
check_prog_verbose minimus2 "minimus2 missing from PATH, download from http://sourceforge.net/apps/mediawiki/amos/index.php?title=AMOS"
check_prog_verbose cd-hit-est "cd-hit-est missing from PATH, required for minimus2 merging, download from http://weizhong-lab.ucsd.edu/cd-hit/"
check_prog_verbose toAmos "toAmos missing from PATH, required for minimus2 merging, download from http://sourceforge.net/apps/mediawiki/amos/index.php?title=AMOS"
check_prog_verbose shuffleSequences_fastq.pl "shuffleSequences_fastq.pl missing from PATH, part of velvet, download from http://www.ebi.ac.uk/~zerbino/velvet"
check_prog_verbose runAssembly "runAssembly missing from PATH, commercial software, required for Newbler merging (NOTE that runAssembly == Newbler)"
check_prog_python_verbose Bio "Bio module not found in python, required for scripts/process-reads/cut-up-fasta.py for Newbler merging, download from http://biopython.org/wiki/Main_Page"
echo "All tests succeeded!!"
