#!/bin/bash -l

## SYNOPSIS ##

# state usage
function usage() {
    echo "Usage: uppmax_raxml.sh -p <phylip> -t <threads> -o <output_dir> -m <standard|bootstrap> [ -a <nt|aa> ]"
    exit
}

# if number or arguments is less than 16, invoke usage function
if (( $# < 8 )); then
    usage;
fi

# defaults
alphabet="aa"

# state options
while getopts ":p:t:a:m:o:" opt; do
    case $opt in
	p) phylip=${OPTARG};;
	t) threads=${OPTARG};;
	a) alphabet=${OPTARG};;
	m) mode=${OPTARG};;
	o) outdir=${OPTARG};;
	*) usage;;
    esac
done

# load tools
module load bioinfo-tools raxml/8.2.4-gcc-mpi

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
if [ "$alphabet" = "nt" ]; then
    model=GTRGAMMA
elif [ "$alphabet" = "aa" ]; then
    model=PROTGAMMALG
fi

echo "Alignment file: $phylip"
echo "Alphabet: $alphabet, Model: $model"
echo "Out directory: $outdir"
echo "Number of threads: $threads"

# do raxml
echo "Running raxml ..."
if [ "$mode" = "standard" ]; then
    raxmlHPC-PTHREADS-AVX \
	-f a -x 12345 -p 12345 -N 100 -m $model \
	-s $phylip -n $runname -w $(pwd)/$outdir/$runname -T $threads &> /dev/null
    echo "Done!"
elif [ "$mode" = "bootstrap" ]; then
    raxmlHPC-PTHREADS-AVX \
	-f d -b 12345 -p 12345 -N 100 -m $model \
	-s $phylip -n $runname -w $(pwd)/$outdir/$runname -T $threads &> /dev/null
    echo "Done!"
fi