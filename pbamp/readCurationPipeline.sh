#!/bin/bash

# Read curation pipeline, as described by Martijn et al, 2019

# usage: readCuratonPipeline.sh <fastq> <sampleName> <oligoFile> <readInformation> <threads>

## oligoFile example:
# primer  CAGCMGCCGCGGTAA CCRAMCTGTCTCACGACG      fwdrev
# primer  CCRAMCTGTCTCACGACG      CAGCMGCCGCGGTAA revfwd
# primer  CAGCMGCCGCGGTAA CAGCMGCCGCGGTAA fwdfwd
# primer  CCRAMCTGTCTCACGACG      CCRAMCTGTCTCACGACG      revrev

## readInformation file example:
# m170201_185651_42203_c101131902550000001823250705311747_s1_p0/54498/ccs 6 0.993807 709
# m170201_185651_42203_c101131902550000001823250705311747_s1_p0/54499/ccs 9 0.99759 1184
# m170201_185651_42203_c101131902550000001823250705311747_s1_p0/54501/ccs 9 0.995953 3960
# m170201_185651_42203_c101131902550000001823250705311747_s1_p0/54502/ccs 4 0.97995 4363
# m170201_185651_42203_c101131902550000001823250705311747_s1_p0/54503/ccs 6 0.992091 3790
# etc ...
## format: <readID>\t<passes>\t<readQuality>\t<readLength>

# Requires that the following softwares are installed:
## seqtk
## python (including Bio, pandas, seaborn, matplotlib)
## mothur
## mafft
## fastx_toolkit
## rnammer
## gnu parallel
## perl (including BioPerl)
## vsearch


# Requires the following custom scripts
## PacBioTrimmer.py
## gff2fasta.pl

# parameters:
fastq_reads=$1
sample=$2
oligos=$3
readinfo=$4
threads=$5


## quality-trim fastq
mkdir 1_fastq2trimfasta

### filter hq reads
echo "Removing reads rq < 0.99 ..."
seqtk subseq \
    $fastq_reads \
    <( awk '$3>0.99 {print $1}' $readinfo ) \
    > 1_fastq2trimfasta/roi.$sample.hq.fastq

### move directory
cd 1_fastq2trimfasta
### remove reads with internal windows of bad read quality
echo "Removing reads with bad quality windows ..."
PacBioTrimmer.py -i roi.$sample.hq.fastq -t 18 -w 30
mv roi.$sample.hq_kept.fastq roi.$sample.rhq.fastq

### fastq.info() -- converts phred 0 to N, converts to FASTA
echo "Running mothur's fastq.info() ..."
mothur --quiet "#fastq.info(fastq=roi.$sample.rhq.fastq, pacbio=T)"

### trim.seqs() -- trims off primers, discards reads with >10homopolymers, length<3000 and >5000, and with 2 or more primer mismatches
echo "Running mothur's trim.seqs() ..."
mothur --quiet "#trim.seqs(fasta=roi.$sample.rhq.fasta, qfile=roi.$sample.rhq.qual, oligos=../oligos.txt, checkorient=T, maxhomop=10, minlength=3000, maxlength=5000, pdiffs=2, keepforward=F, processors=$threads)"

### seqtk subseq -- keep only fwdrev reads
echo "Keeping only fwd-rev reads ..."
seqtk subseq roi.$sample.rhq.trim.fasta <(grep -P "fwdrev|revfwd" roi.$sample.rhq.groups | cut -f 1) > roi.$sample.rhq.trim.fwdrev.fasta

# clean up
echo "Gzip unused files ..."
gzip \
    roi.$sample.hq.fastq \
    roi.$sample.rhq.fasta \
    roi.$sample.rhq.fastq \
    roi.$sample.rhq.groups \
    roi.$sample.rhq.qual \
    roi.$sample.rhq.scrap.fasta \
    roi.$sample.rhq.scrap.qual \
    roi.$sample.rhq.trim.qual \
    roi.$sample.rhq.trim.fasta \
    mothur.*

### move directory
cd ..;
mkdir 2_denovoChimera
cd  2_denovoChimera
ln -s ../1_fastq2trimfasta/roi.$sample.rhq.trim.fwdrev.fasta
### Polarize
echo "Polarizing reads ..."
mafft --adjustdirection --reorder --thread $threads --quiet roi.$sample.rhq.trim.fwdrev.fasta > roi.$sample.rhq.trim.fwdrev.aln
fasta_formatter \
    -i <( sed -r "/^>/! s/\-//g" roi.$sample.rhq.trim.fwdrev.aln ) \
    -o roi.$sample.rhq.trim.fwdrev.pol.fasta \
    -w 60
sed -i -r "/^>/ s/>_R_/>/" roi.$sample.rhq.trim.fwdrev.pol.fasta
### chimera.uchime() -- de novo chimera detection
mothur "#chimera.uchime(fasta=roi.$sample.rhq.trim.fwdrev.pol.fasta, reference=self, chunks=16, abskew=1, chimealns=T)"
mothur "#remove.seqs(fasta=roi.$sample.rhq.trim.fwdrev.pol.fasta, accnos=roi.$sample.rhq.trim.fwdrev.pol.denovo.uchime.accnos)"

### move directory
cd ..;
mkdir 3_rnammer
cd 3_rnammer

ln -s ../2_denovoChimera/roi.$sample.rhq.trim.fwdrev.pol.pick.fasta

# rnammer is really slow, cant multithread, so lets use split parallelize combo
split -d -l 400 roi.$sample.rhq.trim.fwdrev.pol.pick.fasta roi.$sample.rhq.trim.fwdrev.pol.pick.fasta.
# mv each in a new directory, just so that intermediate files don't interfere with each other
for i in roi.$sample.rhq.trim.fwdrev.pol.pick.fasta.*; do
    mkdir $i.dir;
    mv $i $i.dir/;
done

### execute rnammer
parallel -j25 'FOLDER={}; FASTA=${FOLDER%.dir}; cd $FOLDER; rnammer -S bac -multi -gff $FASTA.bac.gff -f $FASTA.bac.fasta $FASTA; cd ..;' ::: roi.$sample.rhq.trim.fwdrev.pol.pick.fasta.*.dir
parallel -j25 'FOLDER={}; FASTA=${FOLDER%.dir}; cd $FOLDER; rnammer -S arc -multi -gff $FASTA.arc.gff -f $FASTA.arc.fasta $FASTA; cd ..;' ::: roi.$sample.rhq.trim.fwdrev.pol.pick.fasta.*.dir

### choose best prediction (arc or bac) based on score (field 6)
cat */*.gff | grep -v "^#" | sort -k1,1 -k9,9 -k6,6 | sort -u -k1,1 -k9,9 > rnammer.$sample.prok.gff
# extract sequences
gff2fasta.pl -g rnammer.$sample.prok.gff -f roi.$sample.rhq.trim.fwdrev.pol.pick.fasta > rnammer.$sample.prok.rRNA.fasta
# glue seqid and desc with _
sed -i '/^>/ s/ /_/' rnammer.$sample.prok.rRNA.fasta
# split 16S and 23S in separate files
seqtk subseq rnammer.$sample.prok.rRNA.fasta <( grep '16s_rRNA' rnammer.$sample.prok.rRNA.fasta | sed "s/>//" ) > rnammer.$sample.prok.16S.fasta
seqtk subseq rnammer.$sample.prok.rRNA.fasta <( grep '23s_rRNA' rnammer.$sample.prok.rRNA.fasta | sed "s/>//" ) > rnammer.$sample.prok.23S.fasta

### select reads that have 16S & 23S predicted
comm -12 <( grep ">" rnammer.$sample.prok.16S.fasta | sed "s/_16s_rRNA//" | sed "s/>//" | sort ) <( grep ">" rnammer.$sample.prok.23S.fasta | sed "s/_23s_rRNA//" | sed "s/>//" | sort) > readsWith16Sand23S.list
seqtk subseq rnammer.$sample.prok.16S.fasta <( cat readsWith16Sand23S.list | while read READ; do grep $READ rnammer.$sample.prok.16S.fasta ; done | sed "s/>//" ) > rnammer.$sample.prok.16S.both.fasta
seqtk subseq rnammer.$sample.prok.23S.fasta <( cat readsWith16Sand23S.list | while read READ; do grep $READ rnammer.$sample.prok.23S.fasta ; done | sed "s/>//" ) > rnammer.$sample.prok.23S.both.fasta

# clean
gzip *.16S.fasta *.23S.fasta *.gff 
rm -r roi.*.dir

### grab ccs surviving ccs reads that have both genes
seqtk subseq ../2_denovoChimera/roi.$sample.rhq.trim.fwdrev.pol.pick.fasta readsWith16Sand23S.list > roi.$sample.rhq.trim.fwdrev.pol.pick.both.fasta

### near-full-length 16S and 23S rRNA gene sequences from the quality-trimmed ccs reads are available in
# rnammer.$sample.prok.16S.both.fasta
# rnammer.$sample.prok.23S.both.fasta

