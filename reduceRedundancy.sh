#!/bin/bash

## this minipipeline aims to reduce the redundancy of a sequence dataset ##
## using UCLUST ##

# usage: #
# $1 = [in -fasta] #
# $2 = [out-dir  ] #
# $3 = [identity ] #

# announce parameters 
echo "FASTA    = $1"
echo "OUTDIR   = $2"
echo "IDENTITY = $3"

# announce start pipeline
echo "####### Clustering ..." $1 "#######"

SEQ=$(basename $1 .fasta)

# sort
uclust \
    --sort $1 \
    --output ${2}/${SEQ}-sort.fasta 

# clust    
uclust \
    --input ${2}/${SEQ}-sort.fasta \
    --uc ${2}/${SEQ}-sort.uc \
    --id $3

# extract seed sequence
uclust \
    --uc2fasta ${2}/${SEQ}-sort.uc \
    --input ${2}/${SEQ}-sort.fasta \
    --output ${2}/${SEQ}-prnd.fasta \
    --types S 

# announce done
echo "####### Done ########"
