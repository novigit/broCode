#!/bin/bash -l

## SYNOPSIS ##
# wrapper bash script to submit spades assemblies to uppmax
# has 4 modes: 
# --only-read-correction (-m hammer) 
# --only_assembler       (-m assembly)
# normal                 (-m standard) (both error correction and assembly)
# --continue             (-m continue)


# state usage
function usage() {
    echo "Usage: uppmax_spades.sh -f <fw.fastq.gz> -r <rv.fastq.gz> -s <singles.fastq.gz> -t <threads> -o <output_dir> -l <memory_limit_in_GB> -k <kmer1,kmer2,kmer3,etc> -m <hammer|assembly|standard|continue>"
    exit
}

# if number or arguments is less than 16, invoke usage function
if [ "$#" -lt "16" ]; then
    usage
fi
echo "the script came here"

# state options
while getopts ":f:r:s:t:o:l:m:k:" opt; do
    case $opt in
	f) fwdReads=${OPTARG};;
	r) rvrReads=${OPTARG};;
	s) unpReads=${OPTARG};;
	t) threads=${OPTARG};;
	o) outpdir=${OPTARG};;
	l) memlimit=${OPTARG};;
	m) mode=${OPTARG};;
	k) kmers=${OPTARG};;
	*) usage;;
    esac
done

# load tools
module load bioinfo-tools spades/3.6.0

# run BayesHammer only
if [ "$mode" = "hammer" ]; then
    spades.py --only-error-correction -1 $fwdReads -2 $rvrReads -s $unpReads -t $threads -o $outpdir -m $memlimit
fi

# run Assembler only
if [ "$mode" = "assembly" ]; then 
    spades.py --sc --only-assembler -1 $fwdReads -2 $rvrReads -s $unpReads -t $threads -o $outpdir -m $memlimit -k $kmers
fi

# run BayesHammer ánd Assembler (standard)
if [ "$mode" = "standard" ]; then
    spades.py --sc -1 $fwdReads -2 $rvrReads -s $unpReads -t $threads -o $outpdir -m $memlimit -k $kmers
fi

# continue run
# if you specify '--continue', the output directory must already exist
if [ "$mode" = "continue" ]; then
    spades.py --continue -o $outpdir
fi