#!/bin/bash

# make taxon to clade mapping file
grep ">" nuclearEncoded29.untreated.aln | sed -e "s/>//" -e "s/$/\tAlpha-nonFE/" > taxa2clade.nu29.map
grep ">" alphamito18.untreated.aln | sed -e "s/>//" -e "s/$/\tAlpha-nonFE/" > taxa2clade.am18.map
grep ">" alphamito24.untreated.aln | sed -e "s/>//" -e "s/$/\tAlpha-nonFE/" > taxa2clade.am24.map
# edit manually, replace Alpha with Alpha-FE, Outgroup or Mito where necessary
cat taxa2clade.am24.map taxa2clade.am18.map taxa2clade.nu29.map | sort -k1,1 -u | sort -k2,2 > taxa2clade.map

# get per taxon per dataset x2scores
# use edited version of script that allows for longer taxon names
for i in *.aln; do perl alignment_pruner.edit.pl --file $i --chi2_test | head -n-2 | awk -v i="${i%.aln}" '{print i,"\t", $2,"\t",$3}' > ${i/aln/x2}; done
cat *.x2 | sed -r "s/ //g" > x2.all.table

# in R do
d<-read.table("x2.all.table", header=F, sep="\t", col.names=c("dataset","taxon","x2score"))

# import clade assignments
map<-read.table("taxa2clade.map", col.names=c("taxa","clade"), header=FALSE)

# merge x2scores with clade info
clades<-merge(d, map, by.x="taxon", by.y="taxa")

# plot all datasets
# set plotting order
clades$dataset<-factor(clades$dataset, levels=c("alphamito24.untreated","alphamito24.f050","alphamito24.stattrim","alphamito18.untreated","alphamito18.f050","nuclearEncoded29.untreated"), ordered=T)
# plot
p<-ggplot(clades, aes(x=dataset, y=x2score)) +
    labs(title="Compositional heterogeneity in different datasets") +
    theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank()) +
    geom_boxplot(outlier.size=NA) +
    geom_jitter(aes(color=clade), position=position_jitter(0.2), size=1) +
    scale_color_manual(values=c("black", "darkgrey","red","blue")) +
    theme(legend.direction = 'horizontal', legend.position ="bottom") +
    geom_signif(test="wilcox.test", test.args=list(alternative="two.sided", paired=TRUE), comparisons=list(c("alphamito24.untreated","alphamito24.f050")), map_signif_level=F, y_position=700) +
    geom_signif(test="wilcox.test", test.args=list(alternative="two.sided", paired=TRUE), comparisons=list(c("alphamito18.untreated","alphamito18.f050")), map_signif_level=F, y_position=700) +
    geom_signif(test="wilcox.test", test.args=list(alternative="two.sided", paired=FALSE), comparisons=list(c("alphamito24.f050","alphamito18.untreated")), map_signif_level=F, y_position=650) +
    geom_signif(test="wilcox.test", test.args=list(alternative="two.sided", paired=FALSE), comparisons=list(c("alphamito24.untreated","alphamito18.untreated")), map_signif_level=F, y_position=750) +
    geom_signif(test="wilcox.test", test.args=list(alternative="two.sided", paired=TRUE), comparisons=list(c("alphamito24.f050","alphamito24.stattrim")), map_signif_level=F, y_position=600) 

ggsave(filename="boxplot.pdf", plot=p, device=cairo_pdf, width=190, height=130, units="mm")
