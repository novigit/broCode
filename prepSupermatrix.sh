#!/bin/bash

# prepares a supermatrix alignment from a bunch of faa's

# state usage
function usage() {
    echo -e "Usage: \n\tprepSupermatrix.sh -f <dir with faas> -o <outdir> -a <linsi|einsi> -p <trimal|bmge> -c <concat_name> -s <separator> \n\n"
    echo "Will replace most X characters with gaps after alignment"
    echo "NOTE: When using special characters like '.' as separators, use two escape charactes. e.g. \\."
    exit
}

# if number or arguments is less than 4, invoke usage function
if [ "$#" -lt "8" ]; then
    usage
    exit
fi

# state options
while getopts ":f:o:a:p:c:s:" opt; do
    case $opt in
	f) faas=${OPTARG};;
	o) outdir=${OPTARG};;
	a) aligner=${OPTARG};;
	p) pruner=${OPTARG};;
	c) concatName=${OPTARG};;
	s) separator=${OPTARG};;
#	t) threads=${OPTARG};;
	*) usage ;;
    esac
done

# prepare out directories
mkdir -p $outdir/1_aln $outdir/2_trim $outdir/3_concat

# align
echo "Aligning clusters in $faas ..."
parallel -j4 "mafft-$aligner --thread 5 --reorder --quiet {} > $outdir/1_aln/{/.}.aln" ::: $faas/*
for i in $outdir/1_aln/*; do sed -i '' -E "/^>/! s/X/-/g" $i; done # replace X characters with gaps


# trim
echo "Trimming clusters in $outdir/1_aln ..."
if [ "$pruner" = "trimal" ]; then
    parallel -j30 "trimal -in {} -out $outdir/2_trim/{/.}.trim.aln -gappyout -fasta" ::: $outdir/1_aln/*
elif [ "$pruner" = "bmge" ]; then
    # # on perun
    # parallel -j30 "java -jar /usr/local/bin/bmge-1.12/BMGE.jar -i {} -of $outdir/2_trim/{/.}.trim.aln -m BLOSUM30 -t AA" ::: $outdir/1_aln/*
    # on macbook
    parallel -j4 "java -Xmx128G -jar /Users/joran/miniconda3/share/bmge-1.12-0/BMGE.jar -i {} -of $outdir/2_trim/{/.}.trim.aln -m BLOSUM30 -t AA" ::: $outdir/1_aln/*
fi


# strip sequence names of everything except species name
for i in $outdir/2_trim/*; do 
    #sed -i '' -E "/^>/ s/$separator.*//" $i;
    sed -i '' -E "/^>/ s/\.\..*//" $i;
done

# concatenate
concatenateRenameAlignment.pl $outdir/2_trim/* > $outdir/3_concat/$concatName.cct.untr.aln
