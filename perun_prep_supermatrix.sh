#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 30
#$ -m bea
#$ -M joran.martijn@dal.ca

# state usage
function usage() {
    echo "Usage: perun_prep_supermatrix.sh -f <fasta_dir> -o <out_dir> -a <linsi|einsi> -p <trimal|bmge> -c <concat_name> -s <separator>"
    echo "Script requires 30 available CPUs"
    exit
}
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

source activate mafft

# align
echo "Aligning clusters in $faas ..."
parallel -j6 "mafft-$aligner --thread 5 --reorder --quiet {} > $outdir/1_aln/{/.}.aln" ::: $faas/*
# replace X characters with gaps
for i in $outdir/1_aln/*; do sed -i -r "/^>/! s/X/-/g" $i; done 

conda deactivate mafft

# trim
echo "Trimming clusters in $outdir/1_aln ..."
if [ "$pruner" = "trimal" ]; then
    parallel -j30 "trimal -in {} -out $outdir/2_trim/{/.}.trim.aln -gappyout -fasta" ::: $outdir/1_aln/*
elif [ "$pruner" = "bmge" ]; then
    source activate bmge
    parallel -j30 "exec java -Xmx128G -jar /scratch2/software/anaconda/envs/bmge/share/bmge-1.12-0/BMGE.jar -i {} -of $outdir/2_trim/{/.}.trim.aln -m BLOSUM30 -t AA" ::: $outdir/1_aln/*
    conda deactivate bmge
fi

# strip sequence names of everything except species name
for i in $outdir/2_trim/*; do 
    sed -i -r "/^>/ s/$separator.*//" $i;
done

# concatenate
export PATH=$PATH:~/repositories/broCode
concatenateRenameAlignment.pl $outdir/2_trim/* > $outdir/3_concat/$concatName.cct.untr.aln