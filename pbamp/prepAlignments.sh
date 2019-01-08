#!/bin/bash

# code for preparing alignments for phylogenies as described by Martijn et al, 2019

# requires the following softwares:
## iqtree-omp
## mafft
## fastx_toolkit
## trimal
## seqtk

### get predicted 16S and 23S sequences from rnammer predictions during read curation pipeline ###
mkdir centroid_rrnas
cd centroid_rrnas/
for sample in PM3 SALA TNS08 P19; do

    # get all predicted 16S and 23S reads
    cp ../../23_reRun_EnvSamples/$sample/3_rnammer/rnammer.$sample.prok.16S.both.fasta . # predicted 16S from all $sample qtrim reads
    cp ../../23_reRun_EnvSamples/$sample/3_rnammer/rnammer.$sample.prok.23S.both.fasta . # predicted 23S from all $sample qtrim reads
    
    # trim '_16s_rRNA' and '_23s_rRNA' from seq headers
    sed -i -r "s/_16s_rRNA//" rnammer.$sample.prok.16S.both.fasta
    sed -i -r "s/_23s_rRNA//" rnammer.$sample.prok.23S.both.fasta
    
    # get centroids
    cp ../../23_reRun_EnvSamples/$sample/4_vsearch/$sample.97p.centroids.fasta . # centroid full $sample qtrim reads

    # select rRNAs from centroids
    seqtk subseq rnammer.$sample.prok.16S.both.fasta <( grep ">" $sample.97p.centroids.fasta | cut -f 3 -d ';' | sed "s/ //g" ) > $sample.centroids.16S.fasta
    seqtk subseq rnammer.$sample.prok.23S.both.fasta <( grep ">" $sample.97p.centroids.fasta | cut -f 3 -d ';' | sed "s/ //g" ) > $sample.centroids.23S.fasta

    # rename centroid sequences, so they include OTU id and OTU size
    grep ">" $sample.97p.centroids.fasta | sed "s/>//" | sed -r "s/([0-9]+)\;size=/\1_s_size=/" | sed -r "s/\;size=.*//" | sed -r "s/\;//g" > $sample.otu2centroid.map
    while read -r otu ccs;
    do
	sed -i -r "/^>/ s|$ccs|${otu}_c_${ccs}|" $sample.centroids.16S.fasta
	sed -i -r "/^>/ s|$ccs|${otu}_c_${ccs}|" $sample.centroids.23S.fasta
    done < $sample.otu2centroid.map
    
    # gzip
    gzip $sample.97p.centroids.fasta rnammer.$sample.prok.16S.both.fasta rnammer.$sample.prok.23S.both.fasta $sample.otu2centroid.map
done
cd ..

### classify reads if bacterial or archaeal ###
mkdir classifySeqs
cd classifySeqs

for sample in PM3 SALA TNS08 P19; do
    cp ../centroid_rrnas/$sample.centroids.16S.fasta .
    cp ../centroid_rrnas/$sample.centroids.23S.fasta .

    # classify
    mothur "#classify.seqs(fasta=$sample.centroids.16S.fasta, reference=silva.nr_v128.align, taxonomy=silva.nr_v128.tax, processors=30)"

    # sort into bacterial and archaeal reads
    seqtk subseq $sample.centroids.16S.fasta <( grep -P "\tBacteria" $sample.centroids.16S.nr_v128.wang.taxonomy | cut -f 1 ) > $sample.centroids.16S.bac.fasta
    seqtk subseq $sample.centroids.16S.fasta <( grep -P "\tArchaea" $sample.centroids.16S.nr_v128.wang.taxonomy | cut -f 1 ) > $sample.centroids.16S.arc.fasta
    seqtk subseq $sample.centroids.23S.fasta <( grep -P "\tBacteria" $sample.centroids.16S.nr_v128.wang.taxonomy | cut -f 1 ) > $sample.centroids.23S.bac.fasta
    seqtk subseq $sample.centroids.23S.fasta <( grep -P "\tArchaea" $sample.centroids.16S.nr_v128.wang.taxonomy | cut -f 1 ) > $sample.centroids.23S.arc.fasta

    # gzip
    gzip $sample.*.accnos $sample.*.taxonomy $sample.*.summary *.logfile $sample.*.16S.fasta $sample.*.23S.fasta
done
cd ..


### 16S+23S phylogeny bacteria ###
mkdir 16s23s_bac
cd 16s23s_bac

# reference 16S bac: 'bac_16S_ref.fasta'
# reference 23S bac: 'bac_23S_ref.fasta'

# pool samples
cat ../classifySeqs/PM3.centroids.16S.bac.fasta ../classifySeqs/SALA.centroids.16S.bac.fasta ../classifySeqs/TNS08.centroids.16S.bac.fasta ../classifySeqs/P19.centroids.16S.bac.fasta > samples.16S.bac.fasta
cat ../classifySeqs/PM3.centroids.23S.bac.fasta ../classifySeqs/SALA.centroids.23S.bac.fasta ../classifySeqs/TNS08.centroids.23S.bac.fasta ../classifySeqs/P19.centroids.23S.bac.fasta > samples.23S.bac.fasta
# check if orientation is correct
# 16S reads should start with smth like 'AGGGT' or 'TACGTAGGGG', in general lots of G's. Not the case! So reverse complement!
# 23S reads should start with smth like 'TCCCTAT'. Not the case! So reverse complement!
seqtk seq -r samples.16S.bac.fasta > samples.16S.bac.rc.fasta # orientation is now OK!
seqtk seq -r samples.23S.bac.fasta > samples.23S.bac.rc.fasta # orientation is now OK!

# align the reference themselves first (may take awhile)
# but these reference alignments can now also be used for 250 bp fragments (see below)
mafft-linsi --thread 30 --maxiterate 1000 --adjustdirection --reorder bac_16S_ref.fasta > bac_16S_ref.aln
mafft-linsi --thread 30 --maxiterate 1000 --adjustdirection --reorder bac_23S_ref.fasta > bac_23S_ref.aln
# check of orientation is correct
## correct orientation 'CAGCMGCCGCGGTAA' (forward primer) and 'CGTCGTGAGACACAGxTxGG' (reverse complement of reverse primer)
# reverse complement is not necessary for 16S!
# reverse complement is not necessary for 23S!

# align 16S and 23S fragments to reference
mafft-linsi --thread 30 --maxiterate 1000 --addfragments samples.16S.bac.rc.fasta bac_16S_ref.aln > ref+samples.16S.bac.aln
mafft-linsi --thread 30 --maxiterate 1000 --addfragments samples.23S.bac.rc.fasta bac_23S_ref.aln > ref+samples.23S.bac.aln

# trim 
# find a gap threshold through trial and error that trims most ambiguous sites while keeping most of the complete reference rRNA genes
# Here: gap threshold 0.3
trimal -in ref+samples.16S.bac.aln -out ref+samples.16S.bac.trim.aln -gt 0.3
trimal -in ref+samples.23S.bac.aln -out ref+samples.23S.bac.trim.aln -gt 0.3

# remove '_R_' from sequence names, introduced by mafft's --adjustdirection
sed -i -r "/^>/ s/>_R_/>/" ref+samples.16S.bac.trim.aln
sed -i -r "/^>/ s/>_R_/>/" ref+samples.23S.bac.trim.aln
# reference has 'U's, reads have 'T's. Make all of them 'T's
sed -i -r "/^>/! s/U/T/gi" ref+samples.16S.bac.trim.aln
sed -i -r "/^>/! s/U/T/gi" ref+samples.23S.bac.trim.aln
# remove centroid name
sed -i -r "/^>/ s/_c_m17.*//" ref+samples.16S.bac.trim.aln
sed -i -r "/^>/ s/_c_m17.*//" ref+samples.23S.bac.trim.aln

# concatenate
concatenateRenameAlignment.pl ref+samples.16S.bac.trim.aln ref+samples.23S.bac.trim.aln > ref+samples.16S23S.bac.trim.cct.aln
# check if orientation alignment is consistent with other alignments for the overall analysis

# gzip
gzip *.fasta *.bac.aln *_ref.aln

### 16S+23S phylogeny archaea ###
cd ..
mkdir 16s23s_arc
cd 16s23s_arc
# reference 16S arc: 'arc_16S_ref.fasta'
# reference 23S arc: 'arc_23S_ref.fasta'

# pool samples
cat ../classifySeqs/PM3.centroids.16S.arc.fasta ../classifySeqs/SALA.centroids.16S.arc.fasta ../classifySeqs/TNS08.centroids.16S.arc.fasta ../classifySeqs/P19.centroids.16S.arc.fasta > samples.16S.arc.fasta
cat ../classifySeqs/PM3.centroids.23S.arc.fasta ../classifySeqs/SALA.centroids.23S.arc.fasta ../classifySeqs/TNS08.centroids.23S.arc.fasta ../classifySeqs/P19.centroids.23S.arc.fasta > samples.23S.arc.fasta
# check if orientation is correct
# 16S reads should start with smth like 'AGGGT' or 'TACGTACCAGCCCC', in general lots of G's and C's. Not the case! So reverse complement!
# 23S reads should start with smth like 'TCCCTAT'. Not the case! So reverse complement!
seqtk seq -r samples.16S.arc.fasta > samples.16S.arc.rc.fasta # orientation is now OK!
seqtk seq -r samples.23S.arc.fasta > samples.23S.arc.rc.fasta # orientation is now OK!

# align the reference themselves first (may take awhile)
# but these reference alignments can now also be used for 250 bp fragments (see below)
mafft-linsi --thread 30 --maxiterate 1000 --adjustdirection --reorder arc_16S_ref.fasta > arc_16S_ref.aln
mafft-linsi --thread 30 --maxiterate 1000 --adjustdirection --reorder arc_23S_ref.fasta > arc_23S_ref.aln
# check of orientation is correct
## correct orientation 'CAGCMGCCGCGGTAA' (forward primer) and 'CGTCGTGAGACACAGxTxGG' (reverse complement of reverse primer)
# reverse complement if necessary
seqtk seq -r arc_16S_ref.aln > arc_16S_ref.rc.aln
seqtk seq -r arc_23S_ref.aln > arc_23S_ref.rc.aln

# # align 16S and 23S fragments to reference
mafft-linsi --thread 30 --maxiterate 1000 --addfragments samples.16S.arc.rc.fasta arc_16S_ref.rc.aln > ref+samples.16S.arc.aln
mafft-linsi --thread 30 --maxiterate 1000 --addfragments samples.23S.arc.rc.fasta arc_23S_ref.rc.aln > ref+samples.23S.arc.aln

# trim 
# find a gap threshold through trial and error that trims most ambiguous sites while keeping most of the complete reference rRNA genes
# Here: gap threshold 0.3
trimal -in ref+samples.16S.arc.aln -out ref+samples.16S.arc.trim.aln -gt 0.3
trimal -in ref+samples.23S.arc.aln -out ref+samples.23S.arc.trim.aln -gt 0.3

# remove '_R_' from sequence names, introduced by mafft's --adjustdirection
sed -i -r "/^>/ s/>_R_/>/" ref+samples.16S.arc.trim.aln
sed -i -r "/^>/ s/>_R_/>/" ref+samples.23S.arc.trim.aln
# reference has 'U's, reads have 'T's. Make all of them 'T's
sed -i -r "/^>/! s/U/T/gi" ref+samples.16S.arc.trim.aln
sed -i -r "/^>/! s/U/T/gi" ref+samples.23S.arc.trim.aln
# remove centroid name
sed -i -r "/^>/ s/_c_m17.*//" ref+samples.16S.arc.trim.aln
sed -i -r "/^>/ s/_c_m17.*//" ref+samples.23S.arc.trim.aln

# concatenate
concatenateRenameAlignment.pl ref+samples.16S.arc.trim.aln ref+samples.23S.arc.trim.aln > ref+samples.16S23S.arc.trim.cct.aln

# gzip
gzip *.fasta *.arc.aln *.bac.aln *.rc.aln *ref.aln

### near-full-length separate 16S and separate 23S bacteria ###
cd ..
mkdir 16s_bac
cd 16s_bac

# get separate 16S and 23S trimmed alignments
ln -s ../16s23s_bac/ref+samples.16S.bac.trim.aln
ln -s ../16s23s_bac/ref+samples.23S.bac.trim.aln
# are in the correct orientation

### near-full-length 16S archaea ###
cd ..
mkdir 16s_arc
cd 16s_arc

# get separate 16S and 23S trimmed alignments
ln -s ../16s23s_arc/ref+samples.16S.arc.trim.aln
ln -s ../16s23s_arc/ref+samples.23S.arc.trim.aln


### 250bp16S bacteria ###
cd ..
mkdir 250bp_16s_bac
cd 250bp_16s_bac

## trim centroids to their first 250 bp
cat ../classifySeqs/PM3.centroids.16S.bac.fasta ../classifySeqs/SALA.centroids.16S.bac.fasta ../classifySeqs/TNS08.centroids.16S.bac.fasta ../classifySeqs/P19.centroids.16S.bac.fasta > samples.16S.bac.fasta # make sure they're all in the right direction

# reverse complement if necessary
seqtk seq -r samples.16S.bac.fasta > samples.16S.bac.rc.fasta

# force upper case notation (fastx toolkit doesn't like lower case)
sed -i -r "/^>/! y/acgt/ACGT/" samples.16S.bac.rc.fasta
# remove ambiguous characters
sed -i -r "/^>/! s/[MRKYWSDVHBXN]//g" samples.16S.bac.rc.fasta
# trim
fastx_trimmer -f 1 -l 250 -i samples.16S.bac.rc.fasta -o samples.16S.bac.rc.250bp.fasta
 
# grab reference alignment
cp ../16s23s_bac/bac_16S_ref.aln.gz . ; gzip -d *.gz;
# reverse complement if necessary
#seqtk seq -r arc_16S_ref.aln > arc_16S_ref.rc.aln
# its in the right direction!
# align 250 bp fragments
mafft-linsi --thread 30 --maxiterate 1000 --addfragments samples.16S.bac.rc.250bp.fasta bac_16S_ref.aln > ref+samples.16S.bac.250bp.aln

# trim
trimal -in ref+samples.16S.bac.250bp.aln -out ref+samples.16S.bac.250bp.trim.aln -gt 0.2

# rename
sed -i -r "/^>/ s/>_R_/>/" ref+samples.16S.bac.250bp.trim.aln 
# replace Us with Ts
sed -i -r "/^>/! s/U/T/gi" ref+samples.16S.bac.250bp.trim.aln 
# remove centroid name
sed -i -r "/^>/ s/_c_m170201.*//" ref+samples.16S.bac.250bp.trim.aln
# gzip
gzip *.fasta *.250bp.aln  *ref.aln

### 250bp16S archaea ###
cd ..
mkdir 250bp_16s_arc
cd 250bp_16s_arc
## trim centroids to their first 250 bp
cat ../classifySeqs/PM3.centroids.16S.arc.fasta ../classifySeqs/SALA.centroids.16S.arc.fasta ../classifySeqs/TNS08.centroids.16S.arc.fasta ../classifySeqs/P19.centroids.16S.arc.fasta > samples.16S.arc.fasta # make sure they're all in the right direction

# reverse complement if necessary
seqtk seq -r samples.16S.arc.fasta > samples.16S.arc.rc.fasta

# force upper case notation (fastx toolkit doesn't like lower case)
sed -i -r "/^>/! y/acgt/ACGT/" samples.16S.arc.rc.fasta
# remove ambiguous characters
sed -i -r "/^>/! s/[MRKYWSDVHBXN]//g" samples.16S.arc.rc.fasta
## == TNS08_o_OTU_91_s_size=1_c_m170201_053342_42203_c101131902550000001823250705311744_s1_p0/125224/ccs == ##
## was found to have too long predicted 16S rRNA genes
## because too long, the first 250 bp were not 16S and resulted in artificially long sequences in phylogenies
## cut 125224 at GGGACG||ATTACTG
## == TNS08_o_OTU_76_s_size=1_c_m170201_053342_42203_c101131902550000001823250705311744_s1_p0/14390/ccs == ##
## was also a long branch in phylogenies
## upon blastn vs nt, it looks like its some sort of artificial sequence, of which several parts mapped to ignisphaera, but also had inserts that did not map to ignisphaera
## remove it
# trim
fastx_trimmer -f 1 -l 250 -i samples.16S.arc.rc.fasta -o samples.16S.arc.rc.250bp.fasta

# grab reference alignment
cp ../16s23s_arc/arc_16S_ref.aln.gz . ; gzip -d *.gz;

# align 250 bp fragments
mafft-linsi --thread 30 --maxiterate 1000 --addfragments samples.16S.arc.rc.250bp.fasta arc_16S_ref.rc.aln > ref+samples.16S.arc.250bp.rc.aln

# trim
trimal -in ref+samples.16S.arc.250bp.rc.aln -out ref+samples.16S.arc.250bp.rc.trim.aln -gt 0.2

# rename
sed -i -r "/^>/ s/>_R_/>/" ref+samples.16S.arc.250bp.rc.trim.aln 
# replace Us with Ts
sed -i -r "/^>/! s/U/T/gi" ref+samples.16S.arc.250bp.rc.trim.aln
# remove centroid name
sed -i -r "/^>/ s/_c_m170201.*//" ref+samples.16S.arc.250bp.rc.trim.aln

# gzip
gzip *.fasta *.rc.aln *_ref.aln
cd ..
# now have the following alignments ready
## ref+samples.16S23S.arc.trim.cct.aln
## ref+samples.16S23S.bac.trim.cct.aln
## ref+samples.16S.arc.trim.aln
## ref+samples.23S.arc.trim.aln
## ref+samples.16S.bac.trim.aln
## ref+samples.23S.bac.trim.aln
## ref+samples.16S.arc.250bp.rc.trim.aln
## ref+samples.16S.bac.250bp.trim.aln
# ready for submitting to iqtree!
