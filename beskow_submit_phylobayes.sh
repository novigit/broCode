#!/bin/bash -l

# documentation
# This script submits the phylobayes with a few standard parameters (running time, number of cores, jobname etc)
# This script calls 'beskow_phylobayes.sh', the one that actually calls phylobayes

# state usage
function usage() {
    echo "Usage: beskow_submit_phylobayes.sh -p <phylip> -m <model> <chains>"
    echo "Will submit each chain, for 23h, with 128 cores each. Note that <chains> has to be the last argument"
    echo "Example: beskow_submit_phylobayes.sh -p concat.phylip -m gtr chain1 chain2 chain3 chain4"
    exit
}

# if number or arguments is less than 5, invoke usage function
if [ "$#" -lt "5" ]; then
    usage
fi

# state options
while getopts ":p:m:" opt; do
    case $opt in
	p) phylip=${OPTARG};;
	m) model=${OPTARG};;
	*) usage;;
    esac
done
shift $(( OPTIND - 1 ))
chains=$*

# submit job
jobname=$(basename $phylip .phylip)
for chain in $chains;
do
    sbatch -A m.2015-1-273 -J $jobname -o $jobname.o -e $jobname.e -t 23:10:00 --nodes=4 --ntasks-per-node=32 beskow_phylobayes.sh $phylip $model $chain 23 128
done