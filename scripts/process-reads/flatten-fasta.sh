#!/bin/bash
set -e
HELPDOC=$( cat <<EOF
Flattens a fasta to have the sequence on one line instead of multiple.

Usage:
    bash [options] inodb-flatten-fasta.sh < fasta file >
Options:
    -h      This help documentation.
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

# Parse arguments
if [ "$#" -ne "1" ]
then
    echo "Invalid number of arguments: 1 needed but $# supplied" >&2
    echo "$HELPDOC"
    exit 1
fi

less $1 | \
awk -v RS='>' \
    -v FS="\n" \
    -v OFS="" \
    -v ORS="" \
    '{ if (NR > 1) { printf ">%s\n",$1; $1=""; printf "%s\n",$0 } }'
