#!/usr/bin/env python

import gffutils
import argparse
import matplotlib.pyplot as plt
import numpy as np

# argument parser
parser = argparse.ArgumentParser(
            description=
        '''
        Parse a GFF3 file, infer intron locations if
        they are not already explicitly defined,
        and calculate relative intron location -relative to the length of the mRNA-
        as well as absolute distances -intron location relative to start of mRNA-

        Then plot histograms for both metrics
        '''
        )
parser.add_argument(
        "-g", "--gff3",
        dest='gff3_file',
        type=str,
        required=True,
        help="Input GFF3 genome file")
parser.add_argument(
        "-t", "--threshold",
        dest='threshold',
        type=int,
        required=False,
        default=1000,
        help="Bin threshold for absolute distance plot")
args = parser.parse_args()

# load gff3 file
db = gffutils.create_db(args.gff3_file, ':memory:', merge_strategy='create_unique')

# infer introns
def transform(f):
    # remove the second exon of the ID
    f.attributes['ID'].pop(1)
    # give the intron feature an intron ID!
    f.attributes['ID'][0] = f.attributes['ID'][0].replace('exon','intron')
    f.source = 'gffutils'
    return f

db.update(db.create_introns(), transform=transform, merge_strategy='create_unique')


# get intron locations
rils = [] # relative
ails = [] # absolute

# iterate over introns
for i in db.features_of_type('intron', order_by=('seqid', 'start')):

    # extract mrna parent
    mrna = db[ i.attributes['Parent'][0] ]
    if i.strand == '+':
        i_rel_start = i.start - mrna.start + 1
    elif i.strand == '-':
        i_rel_start = mrna.end - i.end + 1

    # intron start relative to mrna length
    ril = i_rel_start / len(mrna)
    # intron start relative to mrna start
    ail = i_rel_start
    rils.append(ril)
    ails.append(ail)

    # print(f'{mrna.id}\t{mrna.seqid}\t{mrna.start}\t{ail}')


# set up the 'figure',
## 2 rows, 1 column
fig, axs = plt.subplots(2, 1)

# set up bin boundaries for absolute distance plot
threshold = args.threshold
bin_edges = np.concatenate([ np.arange(0, threshold, 25), [threshold], [max(ails)] ])

# plot absolute distance histogram
axs[0].hist(ails, bins=bin_edges, edgecolor='black')
axs[0].set_xlabel('Intron location relative to mRNA start (bp)')
axs[0].set_xlim(0,threshold+100)
axs[0].set_xticks( np.arange(0, threshold+200, 100) )
xlabels = [int(xtick) if xtick < threshold+100 else f'>{threshold+100}' for xtick in axs[0].get_xticks()]
axs[0].set_xticklabels(xlabels)
axs[0].set_ylabel('Frequency')

# plot relative distance histogram
axs[1].hist(rils, bins=50, edgecolor='black')
axs[1].set_xlabel('Intron location relative to mRNA length')
axs[1].set_ylabel('Frequency')

# adjust space between subplots
plt.subplots_adjust(hspace=1)

# plt.title('Histogram')
plt.savefig('histogram.png')

# show
plt.show()




