#!/bin/bash

# code to evaluate read error rates of reads stemming from the mock community

# usage: assessErrorRates.sh <reads> <reference>
## reads should be supplied in FASTA format

# requires that the following softwares are installed
## blasr
## python3.2
## error_rate_from_sam.py

# parameters:
reads=$1
reference=$2

# map reads to reference
blasr $reads $reference -minMatch 15 -maxMatch 20 -bestn 1 -sam -nproc 10 > ${reads%.fasta}_vs_REF.sam
python3.2 error_rate_from_sam.py -i ${reads%.fasta}_vs_REF.sam -o ${reads%.fasta}_vs_REF.error &> ${reads%.fasta}_vs_REF.error.report
