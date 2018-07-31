#!/bin/bash

# aim: process mock reads again from scratch, but now with the updated pipeline

# parameters:
fastq_reads=$1
sample=$2
oligos=$3
readinfo=$4
threads=$5


mkdir 1_fastq2trimfasta; 
####################
cd 1_fastq2trimfasta
####################

# basename of parameters
fastq=$(basename $fastq_reads)
primers=$(basename $oligos)
info=$(basename $readinfo)

ln -s ../$fastq_reads roi.$sample.fastq
cp ../$oligos .
cp ../$readinfo .


### filter hq reads
echo "Removing reads rq < 0.99 ..."
seqtk subseq \
    roi.$sample.fastq \
    <( awk '$3>0.99 {print $1}' $info ) \
    > roi.$sample.hq.fastq

### remove reads with internal windows of bad read quality
# on Max machine, do; PacBioTrimmer.py -i roi.MOCK.trim.fwdrev.hq.fastq
rsync roi.$sample.hq.fastq $maxlocal; max
cd /local/one/people/jmartijn/
PacBioTrimmer.py -i roi.$sample.hq.fastq
exit
rsync $maxlocal/*$sample*kept.fastq roi.$sample.rhq.fastq

### fastq.info() -- converts phred 0 to N, converts to FASTA
echo "Running mothur's fastq.info() ..."
mothur --quiet "#fastq.info(fastq=roi.$sample.rhq.fastq, pacbio=T)"

### trim.seqs() -- trims off primers, discards reads with >10homopolymers, length<3000 and >5000, and with 2 or more primer mismatches
echo "Running mothur's trim.seqs() ..."
mothur --quiet "#trim.seqs(fasta=roi.$sample.rhq.fasta, qfile=roi.$sample.rhq.qual, oligos=oligos.txt, checkorient=T, maxhomop=10, minlength=3000, maxlength=5000, pdiffs=2, keepforward=F, processors=$threads)"

### seqtk subseq -- keep only fwdrev reads
echo "Keeping only fwd-rev reads ..."
seqtk subseq roi.$sample.rhq.trim.fasta <(grep -P "fwdrev|revfwd" roi.$sample.rhq.groups | cut -f 1) > roi.$sample.rhq.trim.fwdrev.fasta

# clean up
echo "Cleaning up files ..."
rm *.fastq *.groups *.qual *.scrap.fasta *.trim.fasta roi.$sample.rhq.fasta

cd ..;
mkdir 2_denovoChimera
##################
cd 2_denovoChimera
##################

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
mothur "#chimera.uchime(fasta=roi.$sample.rhq.trim.fwdrev.pol.fasta, reference=self, chunks=80, abskew=1, chimealns=T)"
mothur "#remove.seqs(fasta=roi.$sample.rhq.trim.fwdrev.pol.fasta, accnos=roi.$sample.rhq.trim.fwdrev.pol.denovo.uchime.accnos)"

cd ..;
mkdir 3_rnammer
##################
cd 3_rnammer
##################

ln -s ../2_denovoChimera/roi.$sample.rhq.trim.fwdrev.pol.pick.fasta

# rnammer is really slow, cant multithread, so lets use split parallelize combo
split -d -l 400 roi.$sample.rhq.trim.fwdrev.pol.pick.fasta roi.$sample.rhq.trim.fwdrev.pol.pick.fasta.
# mv each in a new directory, just so that intermediate files don't interfere with each other
for i in roi.$sample.rhq.trim.fwdrev.pol.pick.fasta.*; do
    mkdir $i.dir;
    mv $i $i.dir/;
done

# execute rnammer
parallel -j25 'FOLDER={}; FASTA=${FOLDER%.dir}; cd $FOLDER; rnammer -S bac -multi -gff $FASTA.bac.gff -f $FASTA.bac.fasta $FASTA; cd ..;' ::: roi.$sample.rhq.trim.fwdrev.pol.pick.fasta.*.dir
parallel -j25 'FOLDER={}; FASTA=${FOLDER%.dir}; cd $FOLDER; rnammer -S arc -multi -gff $FASTA.arc.gff -f $FASTA.arc.fasta $FASTA; cd ..;' ::: roi.$sample.rhq.trim.fwdrev.pol.pick.fasta.*.dir

# choose best prediction (arc or bac) based on score (field 6)
cat */*.gff | grep -v "^#" | sort -k1,1 -k9,9 -k6,6 | sort -u -k1,1 -k9,9 > rnammer.$sample.prok.gff
# extract sequences
gff2fasta.pl -g rnammer.$sample.prok.gff -f roi.$sample.rhq.trim.fwdrev.pol.pick.fasta > rnammer.$sample.prok.rRNA.fasta
# glue seqid and desc with _
sed -i '/^>/ s/ /_/' rnammer.$sample.prok.rRNA.fasta
# split 16S and 23S in separate files
seqtk subseq rnammer.$sample.prok.rRNA.fasta <( grep '16s_rRNA' rnammer.$sample.prok.rRNA.fasta | sed "s/>//" ) > rnammer.$sample.prok.16S.fasta
seqtk subseq rnammer.$sample.prok.rRNA.fasta <( grep '23s_rRNA' rnammer.$sample.prok.rRNA.fasta | sed "s/>//" ) > rnammer.$sample.prok.23S.fasta

# select reads that have 16S & 23S predicted
comm -12 <( grep ">" rnammer.$sample.prok.16S.fasta | sed "s/_16s_rRNA//" | sed "s/>//" | sort ) <( grep ">" rnammer.$sample.prok.23S.fasta | sed "s/_23s_rRNA//" | sed "s/>//" | sort) > readsWith16Sand23S.list
seqtk subseq rnammer.$sample.prok.16S.fasta <( cat readsWith16Sand23S.list | while read READ; do grep $READ rnammer.$sample.prok.16S.fasta ; done | sed "s/>//" ) > rnammer.$sample.prok.16S.both.fasta
seqtk subseq rnammer.$sample.prok.23S.fasta <( cat readsWith16Sand23S.list | while read READ; do grep $READ rnammer.$sample.prok.23S.fasta ; done | sed "s/>//" ) > rnammer.$sample.prok.23S.both.fasta

cd ..;
mkdir 4_precluster
##################
cd 4_precluster
##################

# grab ccs surviving ccs reads that have both genes
seqtk subseq ../2_denovoChimera/roi.$sample.rhq.trim.fwdrev.pol.pick.fasta ../3_rnammer/readsWith16Sand23S.list > roi.$sample.rhq.trim.fwdrev.pol.pick.both.fasta

# make fasta 99% preClusters
vsearch --cluster_fast roi.$sample.rhq.trim.fwdrev.pol.pick.both.fasta --strand both --id 0.99 --clusters preCluster --threads 20
# add .fasta extension
for i in preCluster*; do mv $i $i.fasta; done
mkdir 1_fasta; mv *[0-9].fasta 1_fasta; 

mkdir 2_mafftQinsi
# # denovo align with mafft
# get list of fasta files with more than 1 sequence
grep -c ">" 1_fasta/*.fasta | grep -vP ':1$' | cut -f 1 -d ':' > files_to_align.txt
parallel -j5 "mafft-qinsi --kimura 1 --thread 5 --quiet {} > 2_mafftQinsi/{/.}.aln" :::: files_to_align.txt

cd ../
mkdir 5_consensus
#################
cd 5_consensus
# ###############

# link alignments (consensus.seqs does not have output option, will write in directory where it found file)
for i in ../4_precluster/2_mafftQinsi/*.aln; do ln -s $i; done

# call consensus for non-singletons
parallel -j30 'aln={}; mothur "#consensus.seqs(fasta=$aln,cutoff=51)"' ::: *.aln

# make sequence names like 'Otu710.16S_size=1' and convert to regular fasta
# non singletons
for i in *.cons.fasta; do 
    otuname=$(basename $i .cons.fasta); 
    otusize=$( grep -c ">" ${i%.cons.fasta}.aln ); 
    echo $i $otuname $otusize; 
    sed -i -r "s/>.*/>${otuname}_i_size=$otusize/" $i;
    fasta_formatter -i <( sed -r "/^>/! s/-//g" $i ) -o ${i%.fasta}.regular.fasta -w 0
done
# remove size=2 clusters
grep size=2$ * | cut -f 1 -d ':' | while read x; do rm $x; done

# cleanup
rm *.aln *cons.fasta *singleton.fasta
mkdir summaries; mv *.summary summaries/
mkdir fastas; mv *.fasta fastas

# pool consensus
cat fastas*/*regular.fasta > $sample.consensus.fasta
# gather reads from singletons and doubletons
grep -c ">" ../4_precluster/1_fasta/* | grep -P "\:1$|\:2$" | cut -f 1 -d ':' | while read x; do cat $x >> $sample.unclustered.fasta; done
cat $sample.*.fasta > $sample.cons+uncl.fasta

cd ..;
mkdir 6_rnammer
##################
cd 6_rnammer
##################

ln -s ../5_consensus/$sample.cons+uncl.fasta
fasta_formatter -i $sample.cons+uncl.fasta -o $sample.cons+uncl.sl.fasta -w 0


# rnammer is really slow, cant multithread, so lets use split parallelize combo
split -d -l 250 $sample.cons+uncl.sl.fasta $sample.cons+uncl.sl.fasta.
# mv each in a new directory, just so that intermediate files don't interfere with each other
for i in $sample.cons+uncl.sl.fasta.*; do
    mkdir $i.dir;
    mv $i $i.dir/;
done

# execute rnammer
parallel -j25 'FOLDER={}; FASTA=${FOLDER%.dir}; cd $FOLDER; rnammer -S bac -multi -gff $FASTA.bac.gff -f $FASTA.bac.fasta $FASTA; cd ..;' ::: $sample.cons+uncl.sl.fasta.*.dir
parallel -j25 'FOLDER={}; FASTA=${FOLDER%.dir}; cd $FOLDER; rnammer -S arc -multi -gff $FASTA.arc.gff -f $FASTA.arc.fasta $FASTA; cd ..;' ::: $sample.cons+uncl.sl.fasta.*.dir

# choose best prediction (arc or bac) based on score (field 6)
cat */*.gff | grep -v "^#" | sort -k1,1 -k9,9 -k6,6 | sort -u -k1,1 -k9,9 > rnammer.$sample.cons+uncl.gff
# extract sequences
gff2fasta.pl -g rnammer.$sample.cons+uncl.gff -f $sample.cons+uncl.sl.fasta > rnammer.$sample.cons+uncl.rRNA.fasta
# glue seqid and desc with _
sed -i '/^>/ s/ /_/' rnammer.$sample.cons+uncl.rRNA.fasta
# split 16S and 23S in separate files
seqtk subseq rnammer.$sample.cons+uncl.rRNA.fasta <( grep '16s_rRNA' rnammer.$sample.cons+uncl.rRNA.fasta | sed "s/>//" ) > rnammer.$sample.cons+uncl.16S.fasta
seqtk subseq rnammer.$sample.cons+uncl.rRNA.fasta <( grep '23s_rRNA' rnammer.$sample.cons+uncl.rRNA.fasta | sed "s/>//" ) > rnammer.$sample.cons+uncl.23S.fasta

cd ..;
mkdir 7_classify
################
cd 7_classify
################

# get silva taxonomy stuff
ln -s ../../4_pipelineTNS08/7_classify/silva.nr_v128.align
ln -s ../../4_pipelineTNS08/7_classify/silva.nr_v128.8mer
ln -s ../../4_pipelineTNS08/7_classify/silva.nr_v128.tax

# classify all otus
ln -s ../6_rnammer/rnammer.TNS08.cons+uncl.16S.fasta
mothur "#classify.seqs(fasta=rnammer.TNS08.cons+uncl.16S.fasta, reference=silva.nr_v128.align, taxonomy=silva.nr_v128.tax, processors=30)"
