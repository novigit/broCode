#!/bin/bash -l

# load parameters
unpReads=$1
threads=$2
outpdir=$3
memlimit=$4

# load tools
module load bioinfo-tools spades/3.5.0

# run BayesHammer
spades.py --only-error-correction -s $unpReads -t $threads -o $outpdir -m $memlimit
