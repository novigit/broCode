#!/bin/bash -l

phylip=$1
model=$2
chain=$3
mode=$4
hours=$5
#profile=$6

# load tools
module load bioinfo-tools gcc/4.6 openmpi/1.4.5 phylobayesmpi/1.8

# start phylobayes chain
if [ "$mode" = "start" ]; then
    mpirun -np 16 pb_mpi -d $phylip -cat -$model $chain &
fi

# continue run
if [ "$mode" = "continue" ]; then
    mpirun -np 16 pb_mpi $chain &
fi

# stop chains before uppmax time runs out
sleep ${hours}h
echo 0 > $chain.run
