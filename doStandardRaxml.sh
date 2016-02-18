#!/bin/bash

# documentation
# Simply does raxml where you dont have to set all the standard settings
# like 100 bootstraps, runname, -x and -p etc

# to do: set a true random number of -p and -x

# state usage
function usage() {
    echo -e "Usage: \n\tdoStandardRaxml.sh -p <phylip> -t <threads> -o <outdir> [ -m <nt|aa> ]\n\n"
    echo "Will replace most illegal characters with underscore in phylip"
    echo "Will use GTRGAMMA if dna, and PROTGAMMALG if protein; Default is PROTGAMMALG"
    echo -e "Final tree files will be in directory with the basename of the phylip\n\n"
    echo "Currently requires <phylip> ends with .phylip. Also number of threads must be > 1"
    exit
}

# if number or arguments is less than 4, invoke usage function
if [ "$#" -lt "6" ]; then
    usage
    exit
fi

# default mode
mode="aa"

# state options
while getopts ":p:t:m:o:" opt; do
    case $opt in
	p) phylip=${OPTARG};;
	t) threads=${OPTARG};;
	m) mode=${OPTARG};;
	o) outdir=${OPTARG};;
	*) usage ;;
    esac
done

# checks
# if file ends with .phylip
# if threads > 1
# if outfile already exists

# prepare outdir
runname=$(basename $phylip .phylip)
mkdir -p $outdir/$runname/

# fix names
sed -i -r "s/[)(:;,]/_/g" $phylip 

# set model
model=PROTGAMMALG # default model
if [ "$mode" = "nt" ]; then
    model=GTRGAMMA
elif [ "$mode" = "aa" ]; then
    model=PROTGAMMALG
fi

echo "Alignment file: $phylip"
echo "Mode: $mode, Model: $model"
echo "Out directory: $outdir"
echo "Number of threads: $threads"

# do raxml
echo "Running raxml ..."
raxmlHPC-PTHREADS-SSE3 \
    -f a -x 12345 -p 12345 -N 100 -m $model \
    -s $phylip -n $runname -w $(pwd)/$outdir/$runname -T $threads &> /dev/null
echo "Done!"
