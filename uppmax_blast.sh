#!/bin/bash -l

# state usage
function usage() {
    echo "Usage: uppmax_blast.sh -q <query> -d <db> -l <taxidlist> -t <threads>"
    exit
}

# load tools
module load bioinfo-tools blast/2.10.1+

# state options
while getopts ":q:d:l:t:" opt; do
    case $opt in
    q) query=${OPTARG};;
    d) db=${OPTARG};;
    l) taxidlist=${OPTARG};;
    t) threads=${OPTARG};;
    *) usage;;
    esac
done

# execute blast
blastp \
    -query $query \
    -db $db \
    -evalue 1e-10 \
    -taxidlist $taxidlist \
    -out ${query}_vs_${db}.blastp \
    -outfmt '6 std qcovs qcovhsp stitle' \
    -max_target_seqs 100000 \
    -num_threads $threads

