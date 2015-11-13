#!/bin/bash


# state usage
function usage() {
    echo "Usage: beskow_submit_phylobayes_convergence.sh -b <burnin> -o <outname> <chains>"
    echo "Example <chains>: 'chain1 chain2 chain3 chain4'. Note that <chains> has to be the last argument"
    echo "<outname> will be the name of the output directory and basename of the output files"
    exit
}

# if number or arguments is less than 16, invoke usage function
if [ "$#" -lt "6" ]; then
    usage
fi

# state options
while getopts ":b:o:" opt; do
    case $opt in
	b) burnin=${OPTARG};;
	o) outname=${OPTARG};;
	*) usage;;
    esac
done

# remove stated options from $*
shift $(( OPTIND - 1 )) # OPTIND's value is the index of the next argument to be processed
chains=$*

# submit job
sbatch -A m.2015-1-273 -J $outname -o $outname.o -e $outname.e -t 3:00 --nodes=1 --ntasks-per-node=1 beskow_phylobayes_convergence.sh $burnin $outname $chains