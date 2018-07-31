#!/usr/local/bin/Rscript

library("ggplot2")
library("reshape2")

args = commandArgs(trailingOnly=TRUE)
matrix<-args[1]

t<-read.table(matrix,header=T)
colnames(t)<-c("Taxon","0","5","10","15","20","25","30","35","40","45","50","55","60","65","70","75","80","85","90","95")
t.long<-melt(t)


pdf("chi2_boxplot.pdf", width=12, height=12)
ggplot(t.long, aes(x=variable,y=value)) + geom_boxplot(alpha=0.75, fatten=3, fill="#BCBCBC") + theme_bw()
dev.off()
