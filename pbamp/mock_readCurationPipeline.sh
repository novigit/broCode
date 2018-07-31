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
echo "Removing reads with bad quality windows ..."
rsync roi.$sample.hq.fastq $maxlocal; max
cd /local/one/people/jmartijn/
PacBioTrimmer.py -i roi.$sample.hq.fastq
exit
rsync $maxlocal/*kept.fastq roi.$sample.rhq.fastq

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
mothur "#chimera.uchime(fasta=roi.$sample.rhq.trim.fwdrev.flip.fasta, reference=self, chunks=80, abskew=1, chimealns=T)"
mothur "#remove.seqs(fasta=roi.$sample.rhq.trim.fwdrev.flip.fasta, accnos=roi.$sample.rhq.trim.fwdrev.flip.denovo.uchime.accnos)"

cd ..;
mkdir 3_rnammer
##################
cd 3_rnammer
##################

ln -s ../2_denovoChimera/roi.$sample.rhq.trim.fwdrev.flip.pick.fasta

# rnammer is really slow, cant multithread, so lets use split parallelize combo
split -d -l 400 roi.$sample.rhq.trim.fwdrev.flip.pick.fasta roi.$sample.rhq.trim.fwdrev.flip.pick.fasta.
# mv each in a new directory, just so that intermediate files don't interfere with each other
for i in roi.$sample.rhq.trim.fwdrev.flip.pick.fasta.*; do
    mkdir $i.dir;
    mv $i $i.dir/;
done

# execute rnammer
parallel -j25 'FOLDER={}; FASTA=${FOLDER%.dir}; cd $FOLDER; rnammer -S bac -multi -gff $FASTA.bac.gff -f $FASTA.bac.fasta $FASTA; cd ..;' ::: roi.$sample.rhq.trim.fwdrev.flip.pick.fasta.*.dir
parallel -j25 'FOLDER={}; FASTA=${FOLDER%.dir}; cd $FOLDER; rnammer -S arc -multi -gff $FASTA.arc.gff -f $FASTA.arc.fasta $FASTA; cd ..;' ::: roi.$sample.rhq.trim.fwdrev.flip.pick.fasta.*.dir

# choose best prediction (arc or bac) based on score (field 6)
cat */*.gff | grep -v "^#" | sort -k1,1 -k9,9 -k6,6 | sort -u -k1,1 -k9,9 > rnammer.$sample.prok.gff
# extract sequences
gff2fasta.pl -g rnammer.$sample.prok.gff -f roi.$sample.rhq.trim.fwdrev.flip.pick.fasta > rnammer.$sample.prok.rRNA.fasta
# glue seqid and desc with _
sed -i '/^>/ s/ /_/' rnammer.$sample.prok.rRNA.fasta
# split 16S and 23S in separate files
seqtk subseq rnammer.$sample.prok.rRNA.fasta <( grep '16s_rRNA' rnammer.$sample.prok.rRNA.fasta | sed "s/>//" ) > rnammer.$sample.prok.16S.fasta
seqtk subseq rnammer.$sample.prok.rRNA.fasta <( grep '23s_rRNA' rnammer.$sample.prok.rRNA.fasta | sed "s/>//" ) > rnammer.$sample.prok.23S.fasta

# select reads that have 16S & 23S predicted
comm -12 <( grep ">" rnammer.MOCK.prok.16S.fasta | sed "s/_16s_rRNA//" | sed "s/>//" | sort ) <( grep ">" rnammer.MOCK.prok.23S.fasta | sed "s/_23s_rRNA//" | sed "s/>//" | sort) > readsWith16Sand23S.list
seqtk subseq rnammer.MOCK.prok.16S.fasta <( cat readsWith16Sand23S.list | while read READ; do grep $READ rnammer.MOCK.prok.16S.fasta ; done | sed "s/>//" ) > rnammer.MOCK.prok.16S.both.fasta
seqtk subseq rnammer.MOCK.prok.23S.fasta <( cat readsWith16Sand23S.list | while read READ; do grep $READ rnammer.MOCK.prok.23S.fasta ; done | sed "s/>//" ) > rnammer.MOCK.prok.23S.both.fasta

# # check error rates at this point
# mkdir 1_error
# cd 1_error
# ln -s ../rnammer.MOCK.prok.16S.both.fasta
# ln -s ../rnammer.MOCK.prok.23S.both.fasta
# cp ../../../9_simplifiedPipelineMock/curatedMockRef/ssuref.mocksubset.addcontaminants.curated.selectpair.fasta .
# cp ../../../9_simplifiedPipelineMock/curatedMockRef/lsuref.mocksubset.addcontaminants.curated.selectpair.fasta .

# mothur "#align.seqs(fasta=rnammer.$sample.prok.16S.both.fasta, reference=ssuref.mocksubset.addcontaminants.curated.selectpair.fasta, flip=T, processors=30)" 
# mothur "#seq.error(fasta=rnammer.$sample.prok.16S.both.align, reference=ssuref.mocksubset.addcontaminants.curated.selectpair.fasta, aligned=T, processors=30)"

# mothur "#align.seqs(fasta=rnammer.$sample.prok.23S.both.fasta, reference=lsuref.mocksubset.addcontaminants.curated.selectpair.fasta, flip=T, processors=30)" 
# mothur "#seq.error(fasta=rnammer.$sample.prok.23S.both.align, reference=lsuref.mocksubset.addcontaminants.curated.selectpair.fasta, aligned=T, processors=30)"


cd ../..;
mkdir 4_precluster
##################
cd 4_precluster
##################

# grab ccs surviving ccs reads that have both genes
seqtk subseq ../2_denovoChimera/roi.MOCK.rhq.trim.fwdrev.flip.pick.fasta ../3_rnammer/readsWith16Sand23S.list > roi.MOCK.rhq.trim.fwdrev.flip.pick.both.fasta

# make fasta 99% preClusters
vsearch --cluster_fast roi.MOCK.rhq.trim.fwdrev.flip.pick.both.fasta --strand both --id 0.99 --clusters preCluster --threads 20 --sizeout
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
# get singleton fastas
comm -3 <( ls ../4_precluster/1_fasta/* ) <( cat ../4_precluster/files_to_align.txt | sed "s/1_fasta/\.\.\/4_precluster\/1_fasta/" ) | while read FILE; do cp $FILE $(basename $FILE .fasta).singleton.fasta; done


# make sequence names like 'Otu710.16S_size=1' and convert to regular fasta
# non singletons
for i in *.cons.fasta; do 
    otuname=$(basename $i .cons.fasta); 
    otusize=$( grep -c ">" ${i%.cons.fasta}.aln ); 
    echo $i $otuname $otusize; 
    sed -i -r "s/>.*/>${otuname}_i_size=$otusize/" $i;
    fasta_formatter -i <( sed -r "/^>/! s/-//g" $i ) -o ${i%.fasta}.regular.fasta -w 0
done
# singletons
for i in *.singleton.fasta; do 
    otuname=$(basename $i .singleton.fasta); 
    echo $i $otuname; 
    sed -i -r "s/>.*/>${otuname}_i_size=1/" $i;
    fasta_formatter -i $i -o ${i%.fasta}.regular.fasta -w 0
done

# cleanup
rm *.aln *cons.fasta *singleton.fasta
mkdir summaries; mv *.summary summaries/
mkdir fastas; mv *.fasta fastas

# pool
cat fastas*/*regular.fasta > $sample.allotus.fasta

cd ..;
mkdir 6_rnammer
##################
cd 6_rnammer
##################

ln -s ../5_consensus/$sample.allotus.fasta

# rnammer is really slow, cant multithread, so lets use split parallelize combo
split -d -l 100 $sample.allotus.fasta $sample.allotus.fasta.
# mv each in a new directory, just so that intermediate files don't interfere with each other
for i in $sample.allotus.fasta.*; do
    mkdir $i.dir;
    mv $i $i.dir/;
done

# execute rnammer
parallel -j25 'FOLDER={}; FASTA=${FOLDER%.dir}; cd $FOLDER; rnammer -S bac -multi -gff $FASTA.bac.gff -f $FASTA.bac.fasta $FASTA; cd ..;' ::: $sample.allotus.fasta.*.dir
parallel -j25 'FOLDER={}; FASTA=${FOLDER%.dir}; cd $FOLDER; rnammer -S arc -multi -gff $FASTA.arc.gff -f $FASTA.arc.fasta $FASTA; cd ..;' ::: $sample.allotus.fasta.*.dir

# choose best prediction (arc or bac) based on score (field 6)
cat */*.gff | grep -v "^#" | sort -k1,1 -k9,9 -k6,6 | sort -u -k1,1 -k9,9 > rnammer.$sample.prok.gff
# extract sequences
gff2fasta.pl -g rnammer.$sample.prok.gff -f $sample.allotus.fasta > rnammer.$sample.prok.rRNA.fasta
# glue seqid and desc with _
sed -i '/^>/ s/ /_/' rnammer.$sample.prok.rRNA.fasta
# split 16S and 23S in separate files
seqtk subseq rnammer.$sample.prok.rRNA.fasta <( grep '16s_rRNA' rnammer.$sample.prok.rRNA.fasta | sed "s/>//" ) > rnammer.$sample.prok.16S.fasta
seqtk subseq rnammer.$sample.prok.rRNA.fasta <( grep '23s_rRNA' rnammer.$sample.prok.rRNA.fasta | sed "s/>//" ) > rnammer.$sample.prok.23S.fasta


# check error rates
cd ..
mkdir 7_checkerrorrates
#######################
cd 7_checkerrorrates
#######################

mkdir 1_rawccs 2_qtrimccs 3_hqconsensus
cp ../../16_fullreadMockRef/SSUitsLSU.mockgenomes.fix.fasta .

cd 1_rawccs/
ln -s ../../0_input/roi.mock.fastq
seqtk seq -a roi.mock.fastq > roi.mock.fasta
blasr roi.mock.fasta ../SSUitsLSU.mockgenomes.fix.fasta -minMatch 15 -maxMatch 20 -bestn 1 -sam -nproc 10 > rawCCS_vs_REF.sam

# make rq vs error plot
cp ../../1_fastq2trimfasta/pb_411_05_ccs_info.txt .
sed -r "s/\/0_[0-9]+\t/\t/" rawCCS_vs_REF.error > rawCCS_vs_REF.error.fix

library(ggplot2)
rq<-read.table("pb_411_05_ccs_info.txt", col.names=c("read","pass","rq","length") )
error<-read.table("rawCCS_vs_REF.error.fix", header=T)
rq2error<-merge(rq,error, by.x="read", by.y="query")
pdf(file="rq_vs_error.pdf")
ggplot(data=rq2error, aes(x=rq, y=error_rate)) + geom_point(alpha=0.1) + labs(x="Read quality", y="Error rate (%)")
dev.off()

cd ..

cd 2_qtrimccs
ln -s ../../4_precluster/roi.MOCK.rhq.trim.fwdrev.flip.pick.both.fasta
blasr roi.MOCK.rhq.trim.fwdrev.flip.pick.both.fasta ../SSUitsLSU.mockgenomes.fix.fasta -minMatch 15 -maxMatch 20 -bestn 1 -sam -nproc 10 > qtrimCCS_vs_REF.sam
cd ..

cd 3_hqconsensus
cp ../../5_consensus/MOCK.allotus.fasta .
seqtk subseq MOCK.allotus.fasta <( grep ">" ../3_consensus/MOCK.allotus.fasta | grep -P -v "_size=1$|size=2$" | sed "s/>//" ) > MOCK.size3orlarger.fasta

mkdir 4_readsNotUsedForConsensus
grep -c ">" ../../4_precluster/1_fasta/* | grep -P ':1$|:2$' | cut -f 1 -d ':' | while read FASTA; do cat $FASTA >> notConsensus.fasta; done

