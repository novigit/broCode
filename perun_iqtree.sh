#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -m bea
#$ -M joran.martijn@dal.ca

# state usage
function usage() {
    echo "Usage: perun_iqtree.sh -s <alignment> -t <threads> [ -p ]"
    echo "-p invokes PMSF"
    exit
}

# state options
while getopts "nfps:t:" opt; do
    case $opt in
	s) alignment=${OPTARG};;
	t) threads=${OPTARG};;
	p) pmsf='triggered';;
	*) usage;;
    esac
done

runname=$(basename $alignment)

# Running PMSF
if [ "$pmsf" = "triggered" ]; then

    echo "PMSF mode was triggered!"

    # create outdirectory
    mkdir $runname.guidetree.out
    mkdir $runname.pmsftree.out

    # generate guidetree
    iqtree-1.6.5 \
        -s $alignment -nt $threads \
        -m LG+G+F -seed 12345 \
        -pre $runname.guidetree.out/$runname.guidetree

    # run PMSF tree (for now fixed to approximation of LG+C60+F+G)
    iqtree-1.6.5 \
        -s $alignment -nt $threads \
        -ft $runname.guidetree.out/$runname.guidetree.treefile \
        -m LG+C60+F+G -b 100 -wbtl -seed 12345 \
        -pre $runname.pmsftree.out/$runname.pmsftree

# Running Normal (default)
else

    echo "Normal mode was triggered!"

    # create outdirectory
    mkdir $runname.out

    iqtree-1.6.5 \
        -s $runname -nt $threads \
        -m TESTNEW -mset LG -madd LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60 \
        -bb 1000 -wbtl -seed 12345 \
        -pre $runname.out/$runname -keep-ident -quiet

fi