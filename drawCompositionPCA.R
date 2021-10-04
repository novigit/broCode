#!/usr/local/bin/Rscript

library("ggplot2")
library("ggfortify")
library("reshape2")
#library("ggbiplot")
library("ggsignif")

args = commandArgs(trailingOnly=TRUE)
amas<-args[1]
mapfile<-args[2]

# read data
data<-read.table(amas,header=TRUE)

# extract amino acid counts
counts<-data[,c(2,seq(6,25,1))]

# calculate total no. of amino acids per taxon
counts$sum<-rowSums(counts[2:21])

# import taxon-to-clade mapping file
map<-read.table(mapfile, col.names=c("taxon","clade"), header=FALSE)

# merge 'counts' with 'map'
clades<-merge(counts, map, by.x="Taxon_name", by.y="taxon")

# convert to frequencies
freqs<-clades[,2:21] / clades$sum


# show contributions of first 5 PC's
pdf(file="pcaVariances.pdf")
pca<-prcomp(freqs)
plot(pca, type="l")
dev.off()

# print to file
pdf(file="compositionPCA.PC1vsPC2.pdf")
# do PCA
autoplot(prcomp(freqs), data=clades, colour="clade", x=1, y=2) + scale_color_manual(values=c("blue","red","orange","grey")) + theme(legend.direction = 'horizontal', legend.position = "bottom")
dev.off()

# print to file
pdf(file="compositionPCA.PC2vsPC3.pdf")
# do PCA
autoplot(prcomp(freqs), data=clades, colour="clade", x=2, y=3) + scale_color_manual(values=c("blue","red","orange","grey")) + theme(legend.direction = 'horizontal', legend.position = "bottom")
dev.off()

# print to file
pdf(file="compositionPCA.PC1vsPC3.pdf")
# do PCA
autoplot(prcomp(freqs), data=clades, colour="clade", x=1, y=3) + scale_color_manual(values=c("blue","red","orange","grey")) + theme(legend.direction = 'horizontal', legend.position = "bottom")
dev.off()

## # print to file
## pdf(file="compositionPCA.PC3vsPC4.pdf")
## # do PCA
## autoplot(prcomp(freqs), data=clades, colour="clade", x=3, y=4) + scale_color_manual(values=c("blue","red","orange","grey")) + theme(legend.direction = 'horizontal', legend.position = "bottom")
## dev.off()

# factor analysis with ggbiplot

# complement freqs table
freqs$taxon<-clades$Taxon_name
freqs$clade<-clades$clade

# pdf(file="biplot.PC1vsPC2.pdf")
# # ggbiplot(pca, choices = c(1,2), groups = freqs$clade) + scale_color_manual(values=c("blue","red","orange","darkgrey")) + theme(legend.direction = 'horizontal', legend.position = "bottom") + geom_text(aes(label=freqs$taxon), size=1)
# ggbiplot(pca, choices = c(1,2), groups = freqs$clade) + scale_color_manual(values=c("blue","red","orange","darkgreen", "purple", "black", "darkgrey")) + theme(legend.direction = 'horizontal', legend.position = "bottom") + geom_text(aes(label=freqs$taxon), size=1)
# dev.off()

# pdf(file="biplot.PC2vsPC3.pdf")
# # ggbiplot(pca, choices = c(2,3), groups = freqs$clade) + scale_color_manual(values=c("blue","red","orange","darkgrey")) + theme(legend.direction = 'horizontal', legend.position = "bottom") + geom_text(aes(label=freqs$taxon), size=1)
# ggbiplot(pca, choices = c(2,3), groups = freqs$clade) + scale_color_manual(values=c("blue","red","orange","darkgreen", "purple", "black", "darkgrey")) + theme(legend.direction = 'horizontal', legend.position = "bottom") + geom_text(aes(label=freqs$taxon), size=1)
# dev.off()

# pdf(file="biplot.PC1vsPC3.pdf")
# # ggbiplot(pca, choices = c(1,3), groups = freqs$clade) + scale_color_manual(values=c("blue","red","orange","darkgrey")) + theme(legend.direction = 'horizontal', legend.position = "bottom") + geom_text(aes(label=freqs$taxon), size=1)
# ggbiplot(pca, choices = c(1,3), groups = freqs$clade) + scale_color_manual(values=c("blue","red","orange","darkgreen", "purple", "black", "darkgrey")) + theme(legend.direction = 'horizontal', legend.position = "bottom") + geom_text(aes(label=freqs$taxon), size=1)
# dev.off()


# make boxplot of distribution of amino acid frequencies
freqs.melt<-melt(freqs)
pdf(file="compositionBoxPlot.pdf")
ggplot(freqs.melt, aes(x=variable, y=value)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color=clade), position=position_jitter(0.2), size=0.75, alpha=0.75) +
  scale_color_manual(values=c("blue","red","orange","darkgrey")) +
  # scale_color_manual(values=c("blue","red","orange","darkgreen", "purple", "black", "darkgrey")) +
  theme(legend.direction = 'horizontal', legend.position = "bottom")
dev.off()

# ############################################
# ### look at particular groups of amino acids (hydrophobic, acidic, polar etc) ###
# ############################################

# clades$hydrophobic<-clades$W + clades$Y + clades$M + clades$I + clades$L + clades$F
# clades$acidic<-clades$D + clades$E
# clades$hydrophilic<-clades$S + clades$N + clades$Q + clades$T

# freqs<-clades[,c(2:21,24:26)] / clades$sum
# freqs$clade<-clades$clade
# freqs$taxon<-clades$Taxon_name

# # violins are set to have the same width, regardless of sample size per group

# ## wilcoxon mann whitney test (non-parametric, not so much assumptions)
# gc<-read.table("GCcontent.csv", header=T)
# freqs.gc<-merge(freqs, gc, by.x="taxon", by.y="Taxon")

# ### GC
# p.gc<-ggplot(data=freqs.gc, aes(x=clade, y=GC)) +
#   geom_boxplot(outlier.size=NA, fill="darkgrey", color="darkred") +
#   geom_jitter(position=position_jitter(0.2), size=0.75, alpha=0.75) +
#   labs(title="%GC") +
#   theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank()) +
#   geom_signif(test="wilcox.test", comparisons = list(c("Methanonatronarchaeia","Haloarchaea")), map_signif_level=F, y_position=74) +
#   geom_signif(test="wilcox.test", comparisons = list(c("Methanonatronarchaeia","MG-IV")), map_signif_level=F, y_position = 50) +
#   geom_signif(test="wilcox.test", comparisons = list(c("MG-IV","rest")), map_signif_level=F, y_position = 22, vjust=2) +
#   geom_signif(test="wilcox.test", comparisons = list(c("Haloarchaea","MG-IV")), map_signif_level=F, y_position = 77) +
#   geom_signif(test="wilcox.test", comparisons = list(c("Haloarchaea","rest")), map_signif_level=F, y_position = 80) +
#   geom_signif(test="wilcox.test", comparisons = list(c("Methanonatronarchaeia","rest")), map_signif_level=F, y_position = 18, vjust=2)
# pdf(file="gcBoxPlot.pdf")
# p.gc
# dev.off()

# ### acidic
# p.acidic<-ggplot(data=freqs, aes(x=clade, y=acidic)) +
#   geom_boxplot(outlier.size=NA, fill="darkgrey", color="darkred") +
#   geom_jitter(position=position_jitter(0.2), size=0.75, alpha=0.75) +
#   labs(title="% Acidic residues (D,E)") +
#   theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank()) +
#   geom_signif(test="wilcox.test", comparisons = list(c("Methanonatronarchaeia","Haloarchaea")), map_signif_level=F, y_position=0.165) +
#   geom_signif(test="wilcox.test", comparisons = list(c("Haloarchaea","rest")), map_signif_level=F, y_position = 0.171) +
#   geom_signif(test="wilcox.test", comparisons = list(c("Haloarchaea","MG-IV")), map_signif_level=F, y_position = 0.168) +
#   geom_signif(test="wilcox.test", comparisons = list(c("Methanonatronarchaeia","MG-IV")), map_signif_level=F, y_position =0.15) +
#   geom_signif(test="wilcox.test", comparisons = list(c("MG-IV","rest")), map_signif_level=F, y_position = 0.087, vjust=2) +
#   geom_signif(test="wilcox.test", comparisons = list(c("Methanonatronarchaeia","rest")), map_signif_level=F, y_position = 0.08, vjust=2)
# pdf(file="acidicBoxPlot.pdf")
# p.acidic
# dev.off()

# ### hydrophobic
# p.hydrophobic<-ggplot(data=freqs, aes(x=clade, y=hydrophobic)) +
#   geom_boxplot(outlier.size=NA, fill="darkgrey", color="darkred") +
#   geom_jitter(position=position_jitter(0.2), size=0.75, alpha=0.75) +
#   labs(title="% Hydrophobic residues (W,Y,M,I,L,F)") +
#   theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank()) +
#   geom_signif(test="wilcox.test", comparisons = list(c("Methanonatronarchaeia","Haloarchaea")), map_signif_level=F, y_position=0.25) +
#   geom_signif(test="wilcox.test", comparisons = list(c("Haloarchaea","MG-IV")), map_signif_level=F, y_position = 0.26) +
#   geom_signif(test="wilcox.test", comparisons = list(c("Haloarchaea","rest")), map_signif_level=F, y_position = 0.27) +
#   geom_signif(test="wilcox.test", comparisons = list(c("Methanonatronarchaeia","rest")), map_signif_level=F, y_position = 0.19, vjust=2) +
#   geom_signif(test="wilcox.test", comparisons = list(c("Methanonatronarchaeia","MG-IV")), map_signif_level=F, y_position =0.24) +
#   geom_signif(test="wilcox.test", comparisons = list(c("MG-IV","rest")), map_signif_level=F, y_position = 0.2, vjust=2)
# pdf(file="hydrophobicBoxPlot.pdf")
# p.hydrophobic
# dev.off()

