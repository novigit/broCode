#!/usr/bin/env python

import argparse
import gffutils

# make command line interface
parser = argparse.ArgumentParser(description="Add intergenic space features to a GFF3 file")
parser.add_argument(
    "-g", 
    dest="gff3_file",
    metavar="GFF3_FILE",
    required=True, 
    help="Input GFF3 file"
)
args = parser.parse_args()

# load GFF3 file into memory
db = gffutils.create_db(args.gff3_file, ':memory:', merge_strategy='create_unique')

# get all gene features
genes = db.features_of_type('gene')
# infer intergenic space features
igss  = db.interfeatures(genes,new_featuretype='intergenic_region')

# transform intergenic space feautures,
# such that there is only one 'ID' attribute
# so from ID=ctg498.gene0459,ctg498.gene046 (2 IDs)
# to      ID=ctg498.gene0459-ctg498.gene046 (1 ID)
# also remove all other attributes for the new igss features if they are there
def transform(f):
    new_id = '-'.join(f.attributes['ID'])
    f.attributes = {}
    f['ID'] = [ new_id ]
    f.source = 'add_igr_space'
    return f

# update db in memory with new intergenic space features
db.update(igss, transform=transform, merge_strategy='create_unique')

# print updated gff record in logical order
for g in db.features_of_type( ('gene','intergenic_region') , order_by=('seqid', 'start')):
    print()
    print(g)
    for f in db.children(g, order_by='start'):
        print(f)

# # print updated gff record in logical order
# for f in db.all_features(order_by=('seqid','start','featuretype'), reverse=True):
#     if f.featuretype == 'gene' or f.featuretype =='intergenic_space':
#         print()
#         print(f)
#     else:
#         print(f)
