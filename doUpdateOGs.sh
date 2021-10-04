#!/bin/bash


function usage() {
    echo "Usage: doUpdateCOGs.sh -i <dir with original faas> -o <dir with updated faas> -f <faas of new genomes> -c <color patterns> -t <threads>"
    exit
}

if [[ $# -lt 8 ]]; then
    usage
fi

# state options
while getopts ":i:o:f:c:t:" opt; do
    case $opt in
	i) original=${OPTARG};;
	o) update=${OPTARG};;
	f) faas=${OPTARG};;
	c) colors=${OPTARG};;
	t) threads=${OPTARG};;
	*) usage ;;
    esac
done

# # align and trim original faas
# mkdir -p $update/1_mafft
# mkdir -p $update/2_trimal
# for originalfaa in $original/*.faa; do

#     faabase=$(basename $originalfaa .faa)

#     echo "Aligning $faabase ..."
#     mafft-linsi \
#     	--thread $threads --quiet $originalfaa \
#     	> $update/1_mafft/${faabase}.aln

#     echo "Trimming alignment ${faabase}.aln ..."
#     trimal \
#     	-in $update/1_mafft/${faabase}.aln \
#     	-out $update/2_trimal/${faabase}_trim.aln \
#     	-gappyout -fasta

# done

# # prepare database of new faa's
# echo "Preparing BLAST database ..."
# cat $faas/*.faa | sed -E "/^>/ s/>/>lcl\|/" > $faas/${faas%\/}.faa.db
# makeblastdb -in $faas/${faas%\/}.faa.db -dbtype prot -parse_seqids &> /dev/null

# do psiblast of originalfaa vs new faa's
mkdir -p $update/3_psiblast
for trimaln in $update/2_trimal/*_trim.aln; do

    faabase=$(basename $trimaln _trim.aln)
    echo "Psiblasting $faabase vs ${faas%\/}..."

    psiblast \
	    -in_msa $trimaln \
	    -db $faas/${faas%\/}.faa.db \
	    -evalue 1e-6 \
	    -outfmt '6 std qcovhsp stitle' \
	    -out $update/3_psiblast/${faabase}_vs_${faas%\/}.psiblast \
	    -num_threads $threads &> /dev/null
done

# retrieve hits and add to originalfaa in updatefaa
mkdir -p $update/4_updated_faas
for psiblast in $update/3_psiblast/*.psiblast; do

    faabase=$(basename $psiblast _vs_${faas%\/}.psiblast)
    echo "Retrieving sequences of psiblast hits of ## $psiblast ## and adding to ## $update/4_updated_faas/$faabase.faa ## ... "
    cp $original/$faabase.faa $update/4_updated_faas/${faabase}_plus_${faas%\/}.faa

    blastdbcmd \
    	-db $faas/${faas%\/}.faa.db \
    	-dbtype prot \
    	-entry_batch <( awk '$3>25 && $13>80 {print $2}' $psiblast | sed -E "s/lcl\|//" | sort -u ) \
    	-outfmt %f \
    	>> $update/4_updated_faas/${faabase}_plus_${faas%\/}.faa
done

# align, trim, fastTree, color
mkdir -p $update/5_mafft $update/6_trimal $update/7_iqtree_fast
for updated_faa in $update/4_updated_faas/*.faa; do

    faabase=$(basename $updated_faa .faa)

    echo -e "\nAligning $faabase.faa ..."
    mafft-linsi --quiet --thread $threads \
    	$updated_faa > $update/5_mafft/${faabase}.aln

    echo "Trimming $faabase.aln ..."
    trimal \
    	-in  $update/5_mafft/${faabase}.aln \
    	-out $update/6_trimal/${faabase}_trim.aln \
    	-gappyout -fasta

    echo "Inferring fast iqtree for $faabase ..."
    iqtree -fast -m LG+G+F -s $update/6_trimal/${faabase}_trim.aln \
    	-pre $update/7_iqtree_fast/${faabase}_trim.aln

    # echo "Coloring $faabase.tree ..."
    # cd $update/7_iqtree_fast/
    # doColorTree.sh -i ${faabase}.tree -p ../../$colors
    # cd ../..

done

echo "Done!"
