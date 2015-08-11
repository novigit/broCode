#!/bin/bash

## SYNOPSIS: ##
# this minipipeline aims to reduce the redundancy of a sequence dataset #
# using UCLUST                                                          #

# state usage
function usage() {
    echo "Usage: reduceRedundancy.sh -i <in.fasta> -o <out.fasta> -p <perc_id>"
    exit
}

if [ $# -lt 6 ]; then
    usage
fi

# state options
while getopts ":i:o:p:" opt; do
    case $opt in
	i) inp=${OPTARG}
	   seq=$(basename $inp .fasta);;
	o) out=${OPTARG};;
	p) idt=${OPTARG};;
	*) usage ;;
    esac
done

# announce start pipeline
echo "####### Clustering ..." $inp "#######"

# sort
uclust \
    --sort $inp \
    --output ${seq}-sort.fasta 

# cluster
uclust \
    --input ${seq}-sort.fasta \
    --uc ${seq}-sort.uc \
    --id $idt \
#    --optimal

# extract seed sequence
uclust \
    --uc2fasta ${seq}-sort.uc \
    --input ${seq}-sort.fasta \
    --output $out \
    --types S 

# edit fasta headers and remove intermediate files
rm ${seq}-sort.uc
rm ${seq}-sort.fasta
sed -i -r "/^>/ s/[0-9]+\|\*\|//" $out

echo "####### Done ########"
