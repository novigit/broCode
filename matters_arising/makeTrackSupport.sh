#!/bin/bash

# do stepwise trimming
## remove spurious 'X' amino acids, inflates chisquare scores
sed -i -r "/^>/! s/X/-/g" lufan_etal_ext_fig_3_alignment.aln
## generate alignments at stepwise trimming steps
parallel -j10 --env _ "perl alignment_pruner.pl --file lufan_etal_ext_fig_3_alignment.aln --chi2_prune f0.{} > lufan_etal_ext_fig_3_alignment.f{}.aln" ::: {000..500..50}

# infer phylogenies under the PMSF approximation of LG+C60+F+G4
for alignment in *.f*.aln; do
    runname=$(basename $alignment)

    # create outdirectory
    mkdir $runname.guidetree.out
    mkdir $runname.pmsftree.out

    # generate guidetree (under LG+G+F)
    iqtree-omp -s $alignment -nt 16 -m LG+G+F -seed 12345 -pre $runname.guidetree.out/$runname.guidetree

    # run PMSF tree 
    iqtree-omp -s $alignment -nt 16 -ft $runname.guidetree.out/$runname.guidetree.treefile -m LG+C60+F+G -b 100 -wbtl -seed 12345 -pre $runname.pmsftree.out/$runname.pmsftree
done

# extract splits (bipartition) counts from IQTREE bootstraps
for i in *.boottrees; do tre_make_splits.pl --file $i --outfile ${i/boottrees/splits}; done

# extract supports for bipartitions of interest
parseSplitcounts.pl -t bipartitionsOfInterest.list -s .aln.pmsftree.splits *.splits | sed "s/lufan_etal_ext_fig_3_alignment\.//g" | cut -f 1-12 > splitsSummary.tsv

# in R do
# store data in 'df'
df<-as.data.frame(t(read.table("splitsSummary.tsv",header=T,row.names=1)))

# make plot
pdf("TrackBPSupports_alphamito_clean.pdf", h=7, w=8)

# leave room for second yaxis
par(mar = c(5,4,4,4) + 0.3) 

# plot RickSister support
plot(df$perc_trim, df$rick_sister, type="b", pch=15, col="black", las=1, xlab="Top biased sites trimmed (%)", ylab="Bootstrap support", ylim=c(0,100))

# add support for other bipartitions
## mito placements 'AlphaSister' and 'AT-rich clade'
points(df$perc_trim, df$alpha_sister, type="b", pch=17, col="darkred")
points(df$perc_trim, df$ATrich_clade, type="b", lty=1, pch=18, col="blue")

## placements of Pelagibacteraceae within 'Core Alphaproteobacteria' and Holobacteraceae within Rhodosprillaceae
points(df$perc_trim, df$Pelagi_Core, type="b", lty=5, pch=2, col="darkgreen")
points(df$perc_trim, df$Holo_Rhodo, type="b", lty=5, pch=2, col="purple")
dev.off()

# add second plot
pdf("TrackX2Score.pdf", h=3, w=8)
plot(df$perc_trim, df$tot_x2score, type="b", bty="n", xlab="", ylab="Total X2-score", col="black")
dev.off()
