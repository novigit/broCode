#!/bin/bash -l
#set -e

## SYNOPSIS ##

# state usage
function usage() {
    echo "Usage: uppmax_iqtree_scratch.sh -s <alignment> -t <threads> [ -p ] [ -f ] [ -n ]"
    exit
}

# if number or arguments is less than 16, invoke usage function
if (( $# < 4 )); then
    usage;
fi

# defaults
#alphabet="aa"

# state options
while getopts "nfps:t:" opt; do
    case $opt in
	s) alignment=${OPTARG};;
	t) threads=${OPTARG};;
	p) pmsf='triggered';;
	f) fast='triggered';;
	n) slow='triggered';;
	*) usage;;
    esac
done

echo "$pmsf" 

# load tools
module load bioinfo-tools iqtree/1.6.5-omp-mpi

# checks
# if file ends with .phylip
# if threads > 1
# if outfile already exists

# fix names
#sed -i -r "s/[)(:;,]/_/g" $alignment

# report stuff
echo "Alignment file: $alignment"
echo "Number of threads: $threads"
echo "Scratch directory: $SNIC_TMP"
if [ "$pmsf" = "triggered" ]; then echo "PMSF mode was triggered!"; fi
if [ "$fast" = "triggered" ]; then echo "Fast IQTREE mode was triggered!"; fi
if [ "$slow" = "triggered" ]; then echo "IQTREE will run with 100 non-parametric bootstraps"; fi


# copy files to $SNIC_TMP, the directory on the scratch disk of the node that the analysis will be run on
picadir=`pwd`
echo "Pica directory: $picadir"

# prepare outdir
runname=$(basename $alignment)
cp    $alignment   $SNIC_TMP
#cp -r $runname.out $SNIC_TMP

# enter $SNIC_TMP
cd $SNIC_TMP

# run iqtree
echo "Running iqtree ..."

if [ "$pmsf" = "triggered" ]; then

    # create outdirectory
    mkdir $runname.guidetree.out
    mkdir $runname.pmsftree.out

    # generate guidetree
    iqtree-omp -s $alignment -nt $threads -m LG+G+F -seed 12345 -pre $runname.guidetree.out/$runname.guidetree

    # run PMSF tree (for now fixed to approximation of LG+C60+F+G)
    iqtree-omp -s $alignment -nt $threads -ft $runname.guidetree.out/$runname.guidetree.treefile -m LG+C60+F+G -b 100 -wbtl -seed 12345 -pre $runname.pmsftree.out/$runname.pmsftree

    # transfer output back to picadir
    cp -r $runname.guidetree.out $picadir
    cp -r $runname.pmsftree.out  $picadir

elif [ "$fast" = "triggered" ]; then

    # create outdirectory
    mkdir $runname.fastiqtree.out

    # run fast iqtree
    iqtree-omp -s $alignment -nt $threads -m TESTNEW -mset LG -madd LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60 -alrt 1000 -wbtl -seed 12345 -pre $runname.fastiqtree.out/$runname -keep-ident -quiet -fast

    # transfer files back to picadir
    cp -r $runname.fastiqtree.out $picadir

elif [ "$slow" = "triggered" ]; then

    # create outdirectory
    #mkdir $runname.npb.nucl.out

    # run slow iqtree
    iqtree-omp -s $alignment -nt $threads -m TESTNEW -mset GTR -b 100 -wbtl -seed 12345 -pre $runname.npb.nucl.out/$runname -keep-ident -quiet

    # transfer files back to picadir
    cp -r $runname.npb.nucl.out $picadir

else

    # create outdirectory
    mkdir $runname.out

    # run iqtree
    echo "iqtree-omp -s $alignment -nt $threads -m TESTNEW -mset LG -madd LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60 -bb 1000 -wbtl -seed 12345 -pre $runname.out/$runname -keep-ident -quiet"
    iqtree-omp -s $alignment -nt $threads -m TESTNEW -mset LG -madd LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60 -bb 1000 -wbtl -seed 12345 -pre $runname.out/$runname -keep-ident -quiet
    
    #transfer files back to picadir
    cp -r $runname.out $picadir

fi

