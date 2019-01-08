#!/bin/bash

# Code for clustering OTUs, as described by Martijn et al, 2019

# usage: clusterOTUs.sh <quality-trimmed reads> <sampleName> <read-quality file> <threads>

# requires:
# vsearch

# parameters:
reads=$1
sample=$2
readinfo=$3

### prep sequence file
mkdir 4_vsearch
sed -r "/^>/ s/\t.*//" 3_rnammer/roi.$sample.rhq.trim.fwdrev.pol.pick.both.fasta > 4_vsearch/$sample.qtrim.fasta
cp 4_vsearch/$sample.qtrim.fasta 4_vsearch/$sample.qtrim.rename.fasta

# format read quality as "size" information in sequence headers
# this lets VSEARCH rank the reads based on read quality prior to clustering
# in the resulting OTUs, the centroid will thus be the one with the highest read quality in the OTU
grep ">" 4_vsearch/$sample.qtrim.fasta | sed "s/>//" | while read READ; do 
    RQ=$(grep $READ $readinfo | cut -f3 -d ' ' | sed "s/0\.//"); 
    echo $READ $RQ; 
    sed -i -r "/^>/ s|>$READ|>${READ};size=$RQ|" 4_vsearch/$sample.qtrim.rename.fasta; 
done
# read qualities don't have trailing 0's. Add them here so conversion to "size" is consistent
sed -i -r "/^>/ s/(size=[0-9]{3})$/\1000/" 4_vsearch/$sample.qtrim.rename.fasta
sed -i -r "/^>/ s/(size=[0-9]{4})$/\100/" 4_vsearch/$sample.qtrim.rename.fasta
sed -i -r "/^>/ s/(size=[0-9]{5})$/\10/" 4_vsearch/$sample.qtrim.rename.fasta

### cluster OTUs
vsearch \
    --cluster_size 4_vsearch/$sample.qtrim.rename.fasta \
    --id 0.97 --strand both --sizeout --sizeorder \
    --relabel ${sample}_o_OTU_ --relabel_keep --threads 20 \
    --centroids 4_vsearch/$sample.97p.centroids.fasta \
    --clusters 4_vsearch/${sample}_o_OTU_

# --cluster_size: order input by size=
# --id 0.97: cluster cutoff at 97% identity
# --strand both: compare both strand directions when making pairwise sequence comparisons
# --sizeout: output size of each 97% OTU in centroid sequence names
# --sizeorder: after clustering, order OTUs by size
# --centroids: output centroids to fasta (here OTU count starts from OTU_1)
# --clusters: output clusters to fasta (starts counting from OTU_0)

### clean dir
gzip 4_vsearch/${sample}_o_OTU_*
