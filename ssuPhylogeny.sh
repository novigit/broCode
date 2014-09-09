#!/bin/bash

## SYNOPSIS: ##
# this minipipeline will infer a maximum likelihood phylogeny for ssu gene #
# using SINA, TRIMAL, RAXML and some optional custom perl scripts          #

# state usage
function usage() {
    echo "Usage: ssuPhylogeny.sh -i [in_fasta] -o [out_tree_name] -d [sina: reference_database] -g [trimal: gap_threshold] -f [figtree_block] -c [color_dictionary] -t [raxml: threads]"
    exit
}

# kill pipeline when there are not sufficient arguments
if [ $# -lt 14 ]; then
    usage
fi

# State options
while getopts ":i:o:d:g:f:c:t:" opt; do
    case $opt in
	i) inp=${OPTARG}
	   seq=$(basename $inp .fasta);;
	o) out=${OPTARG};;
	d) dtb=${OPTARG};;
	g) gpt=${OPTARG};;
	f) fgb=${OPTARG};;
	c) clr=${OPTARG};;
	t) thr=${OPTARG};;
	*) usage ;;
    esac
done


echo -e "###### Start SSU Phylogeny Pipeline for sequences" $inp "######" "\n"

echo -e "\n###### Start SINA alignment ######\n"
# align with SINA
mkdir sina
sina \
    -i $inp \
    -o sina/${seq}.aln \
    --ptdb $dtb \
    --outtype fasta

echo -e "\n###### Start TRIMAL trim ######\n"
# trim with TRIMAL
mkdir trimal
trimal \
    -in sina/${seq}.aln \
    -out trimal/${seq}-trim.phy \
    -gt $gpt \
    -phylip

echo -e "\n###### Start RAXML phylogeny ######\n"
# infer phylogeny with raxml
mkdir raxml
pwd=$(pwd)
raxmlHPC-PTHREADS-SSE3 \
    -f a -x 12345 -p 12345 -N 100 \
    -m GTRGAMMA \
    -s trimal/${seq}-trim.phy \
    -n $out \
    -T $thr \
    -w $pwd/raxml

# convert tree to nexus and color taxa
newick2Nexus.pl -r -i raxml/RAxML_bipartitions.${out} > raxml/${out}.nexus
NexusTreeAddColors.pl \
    -i raxml/${out}.nexus \
    -o raxml/${out}-figrdy.nexus \
    -f $fgb \
    -p $clr
