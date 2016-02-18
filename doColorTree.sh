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
    echo "Usage: doColorTree.sh -i <tree.newick> -p <color_patterns.txt> [ -b <figtree.block> ]"
    echo "If -b is unspecified, will take default figtree.block"
    echo "Output: <.color.nex>, in the same directory as input tree. Open with FigTree"
    exit
}

# if number or arguments is less than 4, invoke usage function
if [ $# -lt 4 ]; then
    usage
fi

# set default figtreeblock
figtreeblock=$(dirname `which $0`)/figtree_block # the same directory as the script

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
upp_convertNewick2Nexus.R $tree

# add figtreeblock and color leafs
NexusTreeAddColors.pl \
    -i $tree.nex \
    -o $tree.color.nex \
    -f $figtreeblock \
    -p $patterns

# remove .nex file (unnecessary)
rm $tree.nex
