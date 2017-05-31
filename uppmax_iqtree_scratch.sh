#!/bin/bash -l
#set -e

## SYNOPSIS ##

# state usage
function usage() {
    echo "Usage: uppmax_iqtree.sh -s <alignment> -t <threads>"
    exit
}

# if number or arguments is less than 16, invoke usage function
if (( $# < 4 )); then
    usage;
fi

# defaults
#alphabet="aa"

# state options
while getopts ":s:t:" opt; do
    case $opt in
	s) alignment=${OPTARG};;
	t) threads=${OPTARG};;
	*) usage;;
    esac
done

# load tools
module load bioinfo-tools iqtree/1.5.3-omp

# checks
# if file ends with .phylip
# if threads > 1
# if outfile already exists

# fix names
sed -i -r "s/[)(:;,]/_/g" $alignment

# report stuff
echo "Alignment file: $alignment"
echo "Number of threads: $threads"
echo "Scratch directory: $SNIC_TMP"

# copy files to $SNIC_TMP, the directory on the scratch disk of the node that the analysis will be run on
picadir=`pwd`
echo "Pica directory: $picadir"

# prepare outdir
runname=$(basename $alignment)
cp    $alignment   $SNIC_TMP
cp -r $runname.out $SNIC_TMP

# enter $SNIC_TMP
cd $SNIC_TMP

# run iqtree
echo "Running iqtree ..."
echo "iqtree-omp -s $alignment -nt $threads -m TESTNEW -mset LG -madd LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60 -bb 1000 -wbtl -seed 12345 -pre $runname.out/$runname -keep-ident -quiet"
iqtree-omp -s $alignment -nt $threads -m TESTNEW -mset LG -madd LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60 -bb 1000 -wbtl -seed 12345 -pre $runname.out/$runname -keep-ident -quiet

# # check if job is properly running
# sleep 1m
# ls -la $SNIC_TMP/$runname.out/              > $picadir/$runname.outputFiles.txt
# cat    $SNIC_TMP/$runname.out/$runname*.log > $picadir/$runname.monitorLog.txt

#mv $runname* $runname.out/

# cp output files back to pica
#sleep 30m
cp -r $runname.out $picadir
