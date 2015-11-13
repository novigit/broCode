#!/bin/bash -l

# set options
phylip=$1
model=$2
chain=$3
hours=$4
threads=$5

# load tools
module load phylobayes/1.5a

# start phylobayes chain
aprun -n $threads pb_mpi -d $phylip -cat -$model -x 10 -1 $chain &

# continue run
# mpirun -np 16 pb_mpi $chain &

# stop chains before time runs out
sleep ${hours}h
echo 0 > $chain.run
