#!/bin/bash

## SYNOPSIS: ##
# this minipipeline aims to reduce the redundancy of a sequence dataset #
# using UCLUST                                                          #

# state usage
function usage() {
    echo "Usage: reduceRedundancy.sh -i [in.fasta] -o [out.fasta] -p [perc_id] -d [out_dir]"
    exit
}

if [ $# -lt 8 ]; then
    usage
fi

# state options
while getopts ":i:o:p:d:" opt; do
    case $opt in
	i) inp=${OPTARG}
	   seq=$(basename $inp .fasta);;
	o) out=${OPTARG};;
	p) idt=${OPTARG};;
	d) dir=${OPTARG};;
	*) usage ;;
    esac
done

# announce start pipeline
echo "####### Clustering ..." $inp "#######"

mkdir $dir
# sort
uclust \
    --sort $inp \
    --output $dir/$seq-sort.fasta 

# cluster
uclust \
    --input $dir/$seq-sort.fasta \
    --uc $dir/$seq-sort.uc \
    --id $idt \
    --optimal \
    --fastapairs $dir/$seq-aln.fasta

# extract seed sequence
uclust \
    --uc2fasta $dir/$seq-sort.uc \
    --input $dir/$seq-sort.fasta \
    --output $dir/$out \
    --types S 

# edit fasta headers
sed -i -r "/^>/ s/[0-9]+\|\*\|//" $dir/$out-prnd.fasta

echo "####### Done ########"
