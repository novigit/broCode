#!/bin/bash -l
#set -e

## SYNOPSIS ##

# state usage
function usage() {
    echo "Usage: uppmax_submit_iqtree_constraint_scratch.sh -s <alignment> -m <model> -c <constraints.list>"
    exit
}

# if number or arguments is less than 16, invoke usage function
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
	c) constraintlist=${OPTARG};;
	p) pmsf='triggered';;
	*) usage;;
    esac
done

# submit constraint searches

if [ "$pmsf" == "triggered" ]; then

    # run constraint search with pmsf model (-p)
    cat $constraintlist | while read topology constraint; do
	echo $constraint > $topology.constraint
	sbatch \
	    -A snic2017-1-592 -J $topology -o $topology.o -e $topology.e -p core -n 8 -t 10-00:00:00 \
	    ~/bin/uppmax_iqtree_constraint_scratch.sh -s $alignment -m $model -c $topology.constraint -p
    done

else

    # run constraint search with lg+c60 model (no -p)
    cat $constraintlist | while read topology constraint; do
    echo $constraint > $topology.constraint
        sbatch \
            -A snic2017-1-592 -J $topology -o $topology.o -e $topology.e -p core -n 8 -t 10-00:00:00 \
            ~/bin/uppmax_iqtree_constraint_scratch.sh -s $alignment -m $model -c $topology.constraint
    done

fi
