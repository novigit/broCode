#!/bin/bash

# state usage
function usage() {
    echo "Usage: doCompleteness.sh -f <fasta> [ -m <bac|arch> ]"
    echo "Default mode: bac"
    echo "Output: .faa, .gff, .hmm, and .report"
    exit
}

# if number or arguments is less than 2, invoke usage function
if [ $# -lt 2 ]; then
    usage
fi

# state options
while getopts ":f:m:" opt; do
    case $opt in
	f) fasta=${OPTARG};;
	m) mode=${OPTARG};;
	*) usage ;;
    esac
done

# set base
name=${fasta%.fasta}

# # do prodigal
# echo "## Predicting ORFs.. ##"
# prodigal-2.60 -q -i $fasta -o $name.gff -a $name.faa

# do mode

# set default
weights=/local/one/tools/micomplete/Bact139.weights
inputhmm=/local/one/tools/micomplete/Bact139.hmm

# if -m is specified
if [ "$mode" = "bac" ]; then
    weights=/local/one/tools/micomplete/Bact139.weights
    inputhmm=/local/one/tools/micomplete/Bact139.hmm
elif [ "$mode" = "arch" ]; then
    weights=/local/one/tools/micomplete/Arch162.weights
    inputhmm=/local/one/tools/micomplete/Arch162.hmm
fi

# do micomplete
echo "## Checking completeness ##"
micomplete.pl \
    -w $weights \
    -h $inputhmm \
    -r $name.completeness.hmm \
    --histogram \
    $fasta \
    &> $name.completeness.report
