#!/bin/bash -l

phylip=$1
model=$2
chain=$3
hours=$4

# load tools
module load bioinfo-tools gcc/4.6 openmpi/1.4.5 phylobayesmpi/1.4f

# start phylobayes chain
# mpirun -np 16 pb_mpi -d $phylip -cat -$model -x 10 -1 $chain &

# continue run
mpirun -np 16 pb_mpi $chain &

# stop chains before uppmax time runs out
sleep ${hours}h
echo 0 > $chain.run
