#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -m bea
#$ -M joran.martijn@dal.ca

alignment=$1
threads=$2

mkdir $alignment.out/
iqtree-1.6.5 -s $alignment -nt $threads -m TESTNEW -mset LG -madd LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60 -bb 1000 -wbtl -seed 12345 -pre $alignment.out/$alignment -keep-ident -quiet
