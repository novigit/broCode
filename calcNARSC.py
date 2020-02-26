#!/usr/bin/env python
import sys
from statistics import mean
from Bio import SeqIO
from Bio.Seq import Seq

# load fasta as file obj
my_fasta = open(sys.argv[1])

# state number of nitrogen
# atoms per amino acid
n_atoms_in_side_chain = {
    'H': 2,
    'K': 1,
    'N': 1,
    'Q': 1,
    'R': 3,
    'W': 1,
}


# define function that calculates narsc per protein
## protein is a Seq() object, 
## mapping a dict with atom counts per aa
def get_narsc_protein(protein, mapping=n_atoms_in_side_chain):
    n_atom_count = 0
    for aa in mapping.keys():
        n_atom_count += protein.count(aa) * mapping.get(aa)
    narsc = n_atom_count / len(protein)
    return narsc

# test function
assert get_narsc_protein(protein=Seq("KNQWRAAAAA"), mapping=n_atoms_in_side_chain)

# parse my_fasta and get list of seq objs
seq_objs = [seq_record.seq for seq_record in SeqIO.parse(my_fasta, format='fasta')]
# get list of narscs and get the mean
print(mean(map(get_narsc_protein, seq_objs)))
# get the mean narsc
# print(mean(all_narscs))
