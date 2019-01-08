#!/bin/bash

# code to generate consensus-ccs sequences, as described by Martijn et al, 2019

# usage: makeConsensusCCS.sh <quality-trimmed ccs reads> <sample> <threads>

## <quality-trimmed ccs reads> are ccs reads that have undergone the following filters
### read-quality > 0.99
### no homopolymers of length 10 or more, between 3000-5000 in length
### are flanked by forward and reverse primers
### are not siamaeras
### are polarized (all reads in the same direction)
### have survived de novo chimera filtering
### encode 16S and 23S genes

# Requires that the following softwares are installed
## vsearch
## mafft
## gnu parallel
## fastx_toolkit

# parameters:
qtrim_reads=$1
sample=$2
threads=$3

### make fasta 99% preClusters
vsearch --cluster_fast $qtrim_reads --strand both --id 0.99 --clusters preCluster --threads $threads --sizeout
# add .fasta extension
for i in preCluster*; do mv $i $i.fasta; done
mkdir 1_fasta; mv *[0-9].fasta 1_fasta; 

### denovo align with mafft
mkdir 2_mafftQinsi
# get list of fasta files with 3 or more sequences
grep -c ">" 1_fasta/*.fasta | grep -vP ':1$|:2$' | cut -f 1 -d ':' > files_to_align.txt
parallel -j5 "mafft-qinsi --kimura 1 --thread 5 --quiet {} > 2_mafftQinsi/{/.}.aln" :::: files_to_align.txt

mkdir 3_consensus
cd 3_consensus

# link alignments (consensus.seqs does not have output option, will write in directory where it found file)
for i in ../2_mafftQinsi/*.aln; do ln -s $i; done

### call consensus
parallel -j30 'aln={}; mothur "#consensus.seqs(fasta=$aln,cutoff=51)"' ::: *.aln

### rename sequences into '>preCluster[0-9]+_i_size=[0-9]+'
for i in *.cons.fasta; do 
    otuname=$(basename $i .cons.fasta); 
    otusize=$( grep -c ">" ${i%.cons.fasta}.aln ); 
    echo $i $otuname $otusize;

    # rename
    sed -i -r "s/>.*/>${otuname}_i_size=$otusize/" $i;

    # remove gaps
    fasta_formatter -i <( sed -r "/^>/! s/-//g" $i ) -o ${i%.fasta}.regular.fasta -w 0
done
# clean up
gzip *.cons.fasta *.summary
cd ..

### add seqs from singleton or doubleton preClusters (no consensus calling)
mkdir 4_singleton_doubleton_seqs
sed -r "s/1_fasta\///g" files_to_align.txt > files_to_exclude.txt # this list contains all 
rsync --exclude-from=files_to_exclude.txt 1_fasta/* 4_singleton_doubleton_seqs
# rename
for i in 4_singleton_doubleton_seqs/*; do
    otuname=$(basename $i .fasta);
    otusize=$( grep -c ">" $i );
    sed -i -r "s/>.*/>${otuname}_i_size=$otusize/" $i;
done

### pool consensus-ccs and non-consensus-ccs seqs
mkdir 5_allSeqs
cat 3_consensus/*.regular.fasta 4_singleton_doubleton_seqs/*.fasta > 5_allSeqs/consensus-ccs_+_non-consensus-ccs_reads.fasta
