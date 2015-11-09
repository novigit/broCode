#!/bin/bash

# state usage
function usage() {
    echo "Usage: doCompleteness.sh -f <fasta> -m <bac|arch>"
    exit
}

# if number or arguments is less than 4, invoke usage function
if [ $# -lt 4 ]; then
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

# do prodigal
echo "## Predicting ORFs.. ##"
prodigal-2.60 -q -i $fasta -o $name.gff -a $name.faa

# do mode
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
    -r $name.hmm \
    --histogram \
    $name.faa \
    &> $name.report
