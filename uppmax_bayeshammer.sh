#!/bin/bash -l

# load parameters
fwdReads=$1
rvrReads=$2
unpReads=$3
threads=$4
outpdir=$5
memlimit=$6

# load tools
module load bioinfo-tools spades/3.5.0

# run BayesHammer
spades.py --only-error-correction -1 $fwdReads -2 $rvrReads -s $unpReads -t $threads -o $outpdir -m $memlimit
