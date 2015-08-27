#!/bin/bash

## SYNOPSIS ##
# takes a tree, converts it to nexus and adds a figtree block and color codes to specific leafs of interest
# output tree is written to same directory as input tree, and gets the extension '.col.nex'

## DEPENDENCIES ##
# assumes that
# 'newick2nexus.pl' and 'NexusTreeAddColors.pl'
# are in your $PATH

# state usage
function usage() {
    echo "Usage: doColorTree.sh -i <tree.newick> -b <figtree.block> -p <color_patterns.txt>"
    exit
}

# if number or arguments is less than 6, invoke usage function
if [ $# -lt 6 ]; then
    usage
fi

# state options
while getopts ":i:b:p:" opt; do
    case $opt in
	i) tree=${OPTARG};;
	b) figtreeblock=${OPTARG};;
	p) patterns=${OPTARG};;
	*) usage ;;
    esac
done

# convert to nexus
newick2nexus.pl -r -i $tree > $tree.nex

# add figtreeblock and color leafs
NexusTreeAddColors.pl \
    -i $tree.nex \
    -o $tree.color.nex \
    -f $figtreeblock \
    -p $patterns

# remove .nex file (unnecessary)
rm $tree.nex
