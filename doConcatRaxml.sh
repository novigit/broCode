#!/bin/bash

## SYNOPSIS ##
# takes a bunch of cogs, aligns them, trims them, makes a concatenation, and infers a raxml tree

## DEPENDENCIES ##
# assumes that 
# 'mafft-linsi, 'trimal', 'concatenateRenameAlignment.pl', and 'raxmlHPC-PTHREADS-SSE3' 
# are in your $PATH

# state usage
function usage() {
    echo "Usage: doConcatRaxml.sh -i <cog_dir> -o <out_dir>"
    exit
}

# if number or arguments is less than 4, invoke usage function
if [ $# -lt 4 ]; then
    usage
fi

# state options
while getopts ":i:o:" opt; do
    case $opt in
	i) faa_dir=${OPTARG};;
	o) out_dir=${OPTARG};;
	*) usage ;;
    esac
done

# prepare outdir
echo -e "###### Prepare output directory $out_dir/ ######"
mkdir -p $out_dir/mafft
mkdir -p $out_dir/trimal
mkdir -p $out_dir/concat
mkdir -p $out_dir/raxml

# do align, trim
for faa in $faa_dir/*.faa;
do
    faa_name=$(basename $faa .faa)
    echo -e "###### Aligning and trimming $faa_name ######"
    mafft-linsi --thread 30 --quiet $faa > $out_dir/mafft/$faa_name.aln
    trimal \
	-in  $out_dir/mafft/$faa_name.aln \
	-out $out_dir/trimal/$faa_name.trim.aln \
	-gt 0.5 \
	-fasta
done

# do concat
echo -e "###### Concatenating ######"
concatenateRenameAlignment.pl $out_dir/trimal/*.trim.aln > $out_dir/concat/concat.trim.aln

# do convert to phylip
trimal \
    -in  $out_dir/concat/concat.trim.aln \
    -out $out_dir/concat/concat.trim.phylip \
    -phylip

# do raxml
echo -e "###### Inferring raxml tree ######"

# do check existence tree with same name
if [ -e $(pwd)/$out_dir/raxml/RAxML_info.concat ];
then
    echo -e "Treefile with the same name already exists!"
fi

raxmlHPC-PTHREADS-SSE3 \
    -f a -x 12345 -p 12345 -N 100 -m PROTGAMMALG \
    -s $out_dir/concat/concat.trim.phylip \
    -n concat -w $(pwd)/$out_dir/raxml -T 30 \
    &> /dev/null
    
echo -e "###### Done! ######"
