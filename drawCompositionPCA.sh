#!/bin/bash

# wrapper script to produce PCA plots for compositional similarities between groups of taxa
# requires R and amas software


# load alignment
alignment=$1
map=$2

# remove 'X' amino acids
sed -i '' -E "/^>/! s/X/-/g" $alignment

# run AMAS
AMAS.py summary -i $alignment -f fasta -d aa -s
summary=${alignment}-seq-summary.txt

# run R script
# script assumes 4 clades/groups of taxa to color
Rscript ~/repositories/broCode/drawCompositionPCA.R $summary $map


#~/scripts/bitbucket/ettemalab/phylogeny/
