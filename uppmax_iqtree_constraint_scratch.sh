#!/bin/bash -l
#set -e

## SYNOPSIS ##

# state usage
function usage() {
    echo "Usage: uppmax_iqtree_constraint_scratch.sh -s <alignment> -m <model> -c <constraint> [ -p ]"
    exit
}

# if number or arguments is less than 8, invoke usage function
if (( $# < 6 )); then
    usage;
fi

# defaults
#alphabet="aa"

# state options
while getopts "ps:m:c:" opt; do
    case $opt in
	s) alignment=${OPTARG};;
	m) model=${OPTARG};;
	c) constraint=${OPTARG};;
	p) pmsf='triggered';;
	*) usage;;
    esac
done

# load tools
module load bioinfo-tools iqtree/1.5.3-omp

# fix names
#sed -i -r "s/[)(:;,]/_/g" $alignment

# report stuff
echo "Alignment file: $alignment"
echo "Model: $model"
echo "Constraint; $constraint"
echo "Scratch directory: $SNIC_TMP"
if [ "$pmsf" = "triggered" ]; then echo "PMSF mode was triggered!"; fi

# copy files to $SNIC_TMP, the directory on the scratch disk of the node that the analysis will be run on
picadir=`pwd`
echo "Pica directory: $picadir"

# prepare outdir
runname=$(basename $alignment .aln)
constraintName=$(basename $constraint .constraint)
cp    $alignment   $SNIC_TMP
cp    $constraint  $SNIC_TMP
#cp -r $runname.out $SNIC_TMP

# enter $SNIC_TMP
cd $SNIC_TMP


if [ "$pmsf" = "triggered" ]; then

    # create outdirs
    mkdir $runname.$constraintName.guidetree.out
    mkdir $runname.$constraintName.pmsftree.out

    # generate guidetree
    iqtree-omp -s $alignment -nt 8 -m LG+G+F -seed 12345 -pre $runname.$constraintName.guidetree.out/$runname.$constraintName.guidetree -keep-ident -quiet

    # run PMSF
    iqtree-omp -s $alignment -nt 8 -g $constraint -ft $runname.$constraintName.guidetree.out/$runname.$constraintName.guidetree.treefile -m $model -seed 12345 -pre $runname.$constraintName.pmsftree.out/$runname.$constraintName.pmsftree -keep-ident -quiet

    # copy results back to picadir
    cp -r $runname.$constraintName.guidetree.out $picadir/
    cp -r $runname.$constraintName.pmsftree.out  $picadir/

else

    #create outdirectory, iqtree doesn't seem to create a new directory when using a not-yet-existing-directory in the -pre flag
    mkdir $runname.$constraintName.out

    #run iqtree
    echo "Running constraint iqtree ..."
    echo "iqtree-omp -s $alignment -nt 8 -m $model -seed 12345 -pre $runname.$constraint.out/$runname.$constraint -keep-ident -quiet"
    iqtree-omp -s $alignment -nt 8 -m $model -g $constraint -seed 12345 -pre $runname.$constraintName.out/$runname.$constraintName -keep-ident -quiet

    #transfer files back to picadir
    cp -r $runname.$constraintName.out $picadir

fi
