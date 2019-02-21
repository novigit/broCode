#!/bin/bash

# Estimate taxonomic compositions of mock and environmental samples, as described by Martijn et al, 2019

# Requires that the following softwares are installed:
## mothur
## R (including ggplot2, RColorBrewer)
## FastX_Toolkit

# parameters:

## run classify.seqs()
### get reference database
ln -s ../24_reRun_Trees/classifySeqs/silva.nr_v128.align
ln -s ../24_reRun_Trees/classifySeqs/silva.nr_v128.8mer
ln -s ../24_reRun_Trees/classifySeqs/silva.nr_v128.tax

## ENVIRONMENTAL SAMPLES
### get all predicted 16S rRNA genes from quality-trimmed reads
for sample in P19 SALA PM3 TNS08; do

    ### link 16S sequences
    ln -s ../23_reRun_EnvSamples/$sample/3_rnammer/rnammer.$sample.prok.16S.both.fasta
    ### also get 250 bp 16S sequences for comparison
    cp ../23_reRun_EnvSamples/$sample/3_rnammer/rnammer.$sample.prok.16S.both.fasta .
    sed -i -r "/^>/! y/acgt/ACGT/" rnammer.$sample.prok.16S.both.fasta
    sed -i -r "/^>/! s/[MRKYWSDVHBXN]//g" rnammer.$sample.prok.16S.both.fasta
    fastx_trimmer -f 1 -l 250 -i rnammer.$sample.prok.16S.both.fasta -o rnammer.$sample.prok.16S.250bp.fasta

    ### run mothur
    mothur "#classify.seqs(fasta=rnammer.$sample.prok.16S.both.fasta, reference=silva.nr_v128.align, taxonomy=silva.nr_v128.tax, relabund=T, processors=30)"
    mothur "#classify.seqs(fasta=rnammer.$sample.prok.16S.250bp.fasta, reference=silva.nr_v128.align, taxonomy=silva.nr_v128.tax, relabund=T, processors=30)"

    ### count number of quality trimmed reads per phylum
    cut -f 2 rnammer.$sample.*.both.*.wang.taxonomy | cut -f 1,2 -d ';' | sed -r "s/\([0-9]+\)//g" | sort | cut -f 2 -d ';' | uniq -c | sed -r "s/^\s+//" | tr ' ' '\t' > $sample.phylumCounts.1000bp.list
    cut -f 2 rnammer.$sample.*.250bp.*.wang.taxonomy | cut -f 1,2 -d ';' | sed -r "s/\([0-9]+\)//g" | sort | cut -f 2 -d ';' | uniq -c | sed -r "s/^\s+//" | tr ' ' '\t' > $sample.phylumCounts.250bp.list

    ### filter out phyla that are less than 0.5% abundant
    readcount=$( grep -c ">" rnammer.$sample.prok.16S.both.fasta )
    cutoffFloat=$( echo "0.005*$readcount" | bc )
    cutoffInt=${cutoffFloat%.*}
    echo $readcount $cutoffFloat $cutoffInt
    awk -v cutoffInt="$cutoffInt" -v sample="${sample}_1000bp" '$1 > cutoffInt {print $0, "\t", sample}' $sample.phylumCounts.1000bp.list > $sample.phylumCounts.1000bp.cutoff.list
    awk -v cutoffInt="$cutoffInt" -v sample="${sample}_250bp" '$1 > cutoffInt {print $0, "\t", sample}' $sample.phylumCounts.250bp.list > $sample.phylumCounts.250bp.cutoff.list

done
cat *.cutoff.list > SAMPLES.phylumCounts.cutoff.list

## visualize composition

# in R do
d<-read.table("SAMPLES.phylumCounts.cutoff.list", header=F, col.names=c("count","phylum","sample"))
library(ggplot2)

# fix levels
archaea<-c("Euryarchaeota","Thaumarchaeota","Aigarchaeota","Crenarchaeota","Bathyarchaeota","Candidate_division_YNPFFA","pMC2A209","Archaea_unclassified")
bacteria<-c("Actinobacteria","Aquificae","Atribacteria","Bacteroidetes","Chloroflexi","Deferribacteres","Gemmatimonadetes","Ignavibacteriae","Nitrospirae","Planctomycetes","Proteobacteria","Spirochaetae","WS2","Chlorobi","Dictyoglomi","Fervidibacteria","Gracilibacteria","Latescibacteria","Omnitrophica","RBG-1_(Zixibacteria)","Bacteria_unclassified")
d$phylum<-factor(d$phylum, levels=c(archaea,bacteria))                                   
                                                                
# fix palette
library(RColorBrewer)
my_palette<-c(brewer.pal(n=7, name="Paired"), "#D3D3D3", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", brewer.pal(n=8, name="Dark2"), brewer.pal(n=9, name="Set1"))

pdf(file="communityProfiles.samples.pdf")
ggplot(d, aes(x=sample, y=count, fill=phylum)) + geom_bar(stat="identity", position="fill") + scale_fill_manual(values=my_palette) + theme(axis.text.x=element_text(angle=90))
dev.off()

## MOCK
ln -s ./23_reRun_EnvSamples/MOCK/3_rnammer/rnammer.MOCK.prok.16S.both.fasta
mothur "#classify.seqs(fasta=rnammer.MOCK.prok.16S.both.fasta, reference=silva.nr_v128.align, taxonomy=silva.nr_v128.tax, relabund=T, processors=30)"
cut -f 2 rnammer.MOCK.*.wang.taxonomy | rev | cut -f 2 -d ';' | rev | sed -r "s/\([0-9]+\)//g" | sort | cut -f 2 -d ';' | uniq -c | sed -r "s/^\s+//" | tr ' ' '\t' | awk '{print $0, "\t", "MOCK"}' > MOCK.genusCounts.list
# in file add species name (we know because we made the mock)
# identified a few extra chimeras! Somehow passed the de novo chimera detection
# they had ambiguous taxonomic assignment, and upon blasting different segments of their sequences against ncbi nt, different mock taxa best hits were found
cp MOCK.genusCounts.list MOCK.genusCounts.noChim.list
# remove them from list
# the read assigned as Pseudomonaceae_unclassified had a good consistent hit vs Pseudomonsa_stutzeri (blast vs nt)
# unclear why this ambiguous assignment
# remove line, and add 1 to Pseudomonas_stutzeri read count

# in R do
d<-read.table("MOCK.genusCounts.noChim.list", header=F, col.names=c("count","genus","species","sample"))
library(ggplot2)

# fix palette
library(RColorBrewer) 
my_palette<-c(brewer.pal(n=12, name="Paired"), brewer.pal(n=8, name="Dark2"), brewer.pal(n=9, name="Set1"), "#7FC97F", "#BEAED4" )
pdf(file="communityProfiles.mock.pdf")
ggplot(d, aes(x=sample, y=count, fill=species)) + geom_bar(stat="identity", position="fill") + scale_fill_manual(values=my_palette)
dev.off()
