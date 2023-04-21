#!/usr/bin/env python

import gffutils
import argparse


# argument parser
parser = argparse.ArgumentParser(
            description=
        '''
        Sort GFF3 file based on contig (seqid), start position and featuretype
        Seperates each gene with an empty newline, so the output GFF3 is
        more easy to read and edit as a human.
        '''
        )
parser.add_argument(
        "-g", "--gff3",
        dest='gff3_file',
        type=str,
        required=True,
        help="Input GFF3 genome file")
args = parser.parse_args()



# load gff file with intergenic regions explicitly defined
db = gffutils.create_db(args.gff3_file, ':memory:', merge_strategy='create_unique')


# print updated gff record in logical order
for g in db.features_of_type('gene', order_by=('seqid', 'start')):
    print()
    print(g)
    for f in db.children(g, order_by='start'):
        print(f)
