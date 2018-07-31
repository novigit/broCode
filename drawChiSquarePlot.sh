#!/bin/bash

# wrapper script to make a chi2plot

# load alignment
alignment=$1
alignment_base=${alignment%.aln}

# remove spurious 'X' amino acids, inflates chisquare scores
sed -i -r "/^>/! s/X/-/g" $alignment

# trim at various steps
parallel -j10 --env _ "perl /home/jmartijn/scripts/bitbucket/ettemalab/phylogeny/alignment_pruner.pl --file $alignment --chi2_prune f0.{} > $alignment_base.f{}.aln" ::: {000..950..50}

# extract per-taxon-chisquare scores
parallel -j10 --env _ "echo {} > $alignment_base.f{}.scores; perl /home/jmartijn/scripts/bitbucket/ettemalab/phylogeny/alignment_pruner.pl --file $alignment_base.f{}.aln --chi2_test | tr -s ' ' | cut -d ' ' -f4 | sed '/^$/d' >> $alignment_base.f{}.scores" ::: {000..950..50}

# build matrix
echo "taxon" > taxa.list
grep ">" $alignment | sed "s/>//" >> taxa.list
paste taxa.list *.scores > $alignment_base.chi2pruning.matrix

# draw plot
Rscript ~/scripts/bitbucket/ettemalab/phylogeny/drawChiSquarePlot.R $alignment_base.chi2pruning.matrix

# clean up
rm *.scores $alignment_base.f*.aln taxa.list $alignment_base.chi2pruning.matrix
