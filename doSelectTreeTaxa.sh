#!/bin/bash

## SYNOPSIS ##
# takes a figtree annotated nexus file, and returns a list of taxa that have been colored in figtree
# thus the actual selection of taxa occurs within FigTree

# state usage
function usage() {
    echo "Usage: doSelectTreeTaxa.sh -i <tree.nexus> -o <selected_taxa.list>"
    exit
}

# if number or arguments is less than 4, invoke usage function
if [ $# -lt 4 ]; then
    usage
fi

# state options
while getopts ":i:o:" opt; do
    case $opt in
	i) tree=${OPTARG};;
	o) list=${OPTARG};;
	*) usage ;;
    esac
done

# make the list
grep "color=" $tree | sed -r "s/\[.*$//" | sed -r "s/\s+//g" > $list
