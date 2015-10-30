#!/bin/bash -l

## set parameters
fw=$1
rv=$2
il_fasta=$3
outdir=$4
cores=$5

## load tools
module load bioinfo-tools IDBA/1.1.1-384

## execute assembly
# convert to fasta and make interleaved
fq2fa --merge <(zcat $fw) <(zcat $rv) $il_fasta

# run assembly
idba_ud --num_threads $cores -r $il_fasta -o $outdir
