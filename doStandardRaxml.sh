#!/bin/bash

# documentation
# Simply does raxml where you dont have to set all the standard settings
# like 100 bootstraps, runname, -x and -p etc

# to do: set a true random number of -p and -x

# state usage
function usage() {
    echo "Usage: doStandardRaxml.sh -p <phylip> -t <threads> [ -m <nt|aa> ]"
    echo "Will use GTRGAMMA if dna, and PROTGAMMALG if protein; Default is PROTGAMMALG"
    echo "Final tree files will be in directory with the basename of the phylip"
    exit
}

# if number or arguments is less than 4, invoke usage function
if [ "$#" -lt "4" ]; then
    usage
    exit
fi

# state options
while getopts ":p:t:m:" opt; do
    case $opt in
	p) phylip=${OPTARG};;
	t) threads=${OPTARG};;
	m) mode=${OPTARG};;
	*) usage ;;
    esac
done

# prepare outdir
runname=$(basename $phylip .phylip)
outdir=$(pwd)/$runname
mkdir $outdir

# set model
model=PROTGAMMALG # default model
if [ "$mode" = "nt" ]; then
    model=GTRGAMMA
elif [ "$mode" = "aa" ]; then
    model=PROTGAMMALG
fi

# do raxml
raxmlHPC-PTHREADS-SSE3 \
    -f a -x 12345 -p 12345 -N 100 -m $model \
    -s $phylip -n $runname -w $outdir -T $threads \
    &> /dev/null
