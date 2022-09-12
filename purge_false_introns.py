#!/usr/bin/env python

import gffutils
import pysam
import orfipy_core
import re
from Bio import SeqIO
from statistics import mean

# load gff file into memory
db = gffutils.create_db(
    data='test.simple.gff3',
    # dbfn='test.simple.gff3.db',
    dbfn=':memory:',
    merge_strategy="error"
)

######################## - ADD INTRON FEATURES - ###############################

# infer and add intron features
introns = list(db.create_introns())

# edit ID attributes of introns
## intron features have ID=exon1,exon2 as format
## this breakse .update(), so we need to adjust their IDs
for i in introns:

    # remove the second exon of the ID
    i.attributes['ID'].pop(1)

    # give the intron feature an intron ID!
    i.attributes['ID'][0] = i.attributes['ID'][0].replace('exon','intron')

    # change source field from 'gffutils_derived'
    # to 'gffutils'
    i.source = 'gffutils'

# update db with new introns
db.update(introns)

######################## - IDENTIFY GENES WITH FALSE INTRONS - #################

# load bam file
bamfile = pysam.AlignmentFile('rnaseq_vs_masked_ergo_cyp_genome.sort.bam', mode='rb')

# initiate list of genes with false introns
genes_with_false_introns = []

# loop over introns, check which are real or false, and update intron features
updated_introns = []
introns = db.features_of_type(featuretype='intron', limit='tig00000012:1-400000')
for i in introns:
    
    # start list of spliced / all aligned reads ratio's
    ratios = []

    # i.start and i.end are 1-indexed
    # pysam works with 0-indexing
    # so to get coverage for position x (1-indexed),
    # we need to ask for position x-1 (0-indexed)
    for col in bamfile.pileup(contig=i.seqid, start=i.start-1, stop=i.end, truncate=True):

        # skip col if no reads are aligned
        if col.get_num_aligned() == 0:
            continue
        
        # start spliced alignment counter
        spliced = 0

        # iterate over aligned reads at col
        for read in col.pileups:
            # add 1 to spliced alignment counter if read has 'N' here
            # in the CIGAR string. This indicates the read is spliced here
            if read.is_refskip:
                spliced = spliced + 1

        # calculate ratio (spliced / all aligned reads)
        ratio = spliced / col.get_num_aligned()

        # append ratio to list of ratios
        ratios.append(ratio)
            
    # if no reads map to this intron at all,
    # we can not be sure it is not an intron
    # so we trust the prediction software and keep it
    if len(ratios) == 0:
        ratios = [1.0]

    # if ratio < 0.80,
    # intron is most likely false
    if mean(ratios) < 0.50:

        # mark intron as false
        i.attributes['call'] = 'false'

        # add gene to list of 
        # genes with false introns
        parent_gene = db.parents(i, featuretype='gene')
        for p in parent_gene:
            if p not in genes_with_false_introns:
                genes_with_false_introns.append(p)


    # if ratio > 0.80,
    # intron is most likely true
    elif mean(ratios) >= 0.50:
        
        # mark intron as true
        i.attributes['call'] = 'true'

    # add to updated_introns
    updated_introns.append(i)

    # delete original intron
    db.delete(i)

# update db with updated intron features
db.update(updated_introns)


######################## - UPDATE FEATURES OF GENES WITH FALSE INTRONS - ##

# load FASTA file
fasta = SeqIO.index('ergo_cyp_genome.fasta','fasta')

new_features = []
possible_phases = [0,2,1]

# generate new feature objects from template
def new_feature_from_template(gene, ftype, strand):
    if ftype == 'gene':
        gene_id = gene.id
        attrs = gene.attributes
    else:
        gene_id = gene.id.replace('gene','mRNA')
        attrs = 'ID='+gene_id+'.'+ftype+'01;Parent='+gene_id
    # define template feature
    feature_template = gffutils.feature.Feature(
            seqid=gene.seqid, 
            source='gffutils', 
            featuretype=ftype, 
            start=gene.start, 
            end=gene.end, 
            score='.',
            strand=strand,  
            frame='.',
            attributes=attrs
    )
    return feature_template


# given a sequence and strand,
# return a list of non-overlapping orfs
def get_non_overlapping_orfs(seq, strand):
    # convert strand symbol
    if strand == '+':
        strand = 'f'
    else:
        strand = 'r'
    # find orfs in seq
    ## orfs() returns a list of tuples
    all_orfs = orfipy_core.orfs(
        seq,
        minlen=300,
        maxlen=10000, 
        starts=['ATG'], 
        strand=strand, 
        include_stop=True
    )
    # first sort orfs by length
    all_orfs.sort(key=lambda x : x[1]-x[0], reverse=True) 
    # then select non-overlapping orfs,
    # preferring the largest ones
    occupied_regions = []
    selected_orfs = []
    for orf in all_orfs:
        # get coordinates, frame and length of orf
        start = orf[0] 
        stop  = orf[1]
        # always select the first orf, 
        # i.e. the largest orf
        if len(selected_orfs) == 0:
            occupied_regions.append([start,stop])
            selected_orfs.append(orf)
            continue
        # check if current orf overlaps with
        # any of the already selected orfs
        for region in occupied_regions:
            # move to next orf if any region overlaps 
            if max(region[0], start) < min(region[1], stop):
                break
        # this else block is only executed if the above
        # for loop is exited normally (i.e. doesn't break)
        else:
            # if it survived the overlap checks, select orf
            occupied_regions.append([start,stop])
            selected_orfs.append(orf)
    return selected_orfs
    



# loop over genes with false introns
for g in genes_with_false_introns:

    # define region and strand on which we have to operate
    ## get start coordinate
    region_start = g.start
    ## get end coordinate (i.e. start of next gene - 1)
    gene_number = int( g.id.split('gene')[1] )
    gene_contig = str( g.id.split('gene')[0] )
    next_gene_id = gene_contig + "gene%04d" % (gene_number + 1)
    region_end = db[next_gene_id].start - 1
    ## extract region sequence
    ### note that slicing is 0-indexed, and gff coordinates 1-indexed
    contig = fasta[g.seqid].seq
    region_seq = str( contig[region_start-1:region_end] )
    ## get strand
    strand = g.strand

    # determine if there are any true introns in this gene
    ## get true intron borders
    introns = list( db.children(g, featuretype='intron', order_by=('seqid','start','attributes')) )
    true_introns = filter(lambda x : x.attributes['call'][0] == 'true', introns)
    intron_borders = [ x for t in true_introns for x in [t.start, t.end] ]
    ## delete false introns
    false_introns = filter(lambda x : x.attributes['call'][0] == 'false', introns)
    db.delete(false_introns)

    # create new template gene and mRNA feature
    new_gene = new_feature_from_template(g, 'gene', strand)
    new_mRNA = new_feature_from_template(g, 'mRNA', strand)

    # apply naive orf finder if there are no introns
    if len(intron_borders) == 0:
        print(g.id)
        new_orfs = get_non_overlapping_orfs(region_seq, strand)
        i=1
        for orf in new_orfs:
            new_exon = new_feature_from_template(g, 'exon', strand)
            # adjust coordinates
            new_exon.start = region_start + orf[0]
            new_exon.end   = region_start + orf[1] - 1
            # update name
            new_exon.attributes['ID'][0] = new_exon.attributes['ID'][0] + '.' + str(i)
            new_exon.id = new_exon.attributes['ID'][0]
            print(new_exon)
            i = i+1


        # # traverse gene sequence in steps of 3
        # # or by jumping over true introns until 
        # # start of next gene and create exons along the way
        # n = g.start - 1         # g.start is 1-indexed, n should be 0-indexed
        # exon_start = g.start    # initialize exon_start, this will be updated as we traverse the seq
        # exon_number = 1         # start exon counter
        # while n < next_gene_start - 1:

        #     # if we encounter a true intron
        #     if len(intron_borders) != 0 and n >= intron_borders[0]:
        #         # create new exon and CDS features
        #         new_exon = new_feature_from_template(g, 'exon')
        #         new_CDS  = new_feature_from_template(g, 'CDS')

        #         # update new_exon properties
        #         new_exon.start = exon_start
        #         new_exon.end   = intron_borders[0] - 1      # intron_borders[0] is the start of the true intron
        #         pattern        = "exon%02d" % exon_number 
        #         new_exon.attributes['ID'][0] = new_exon.attributes['ID'][0].replace('exon01', pattern)
        #         new_exon.id = new_exon.attributes['ID'][0]

        #         # update new_CDS properties
        #         new_CDS.start = exon_start
        #         new_CDS.end   = intron_borders[0] - 1      # intron_borders[0] is the start of the true intron
        #         new_CDS.frame = str(phase)
        #         pattern       = "CDS%02d" % exon_number 
                # new_CDS.attributes['ID'][0] = new_CDS.attributes['ID'][0].replace('CDS01', pattern)
                # new_CDS.id = new_CDS.attributes['ID'][0]

                # # update exon counter
                # exon_number += 1

                # # add new exon to new_features
                # new_features.extend([new_exon, new_CDS])

                # # jump over intron to start of next exon and adjust phase
                # codon_remainder = (new_exon.end - new_exon.start + 1) % 3
                # phase = possible_phases[codon_remainder]
                # exon_start = intron_borders[1] + 1          # intron_borders[1] is the end of the true intron
                # n = exon_start + phase - 1

                # # delete first pair of start and end of true introns from intron_borders
                # # to ensure in next loop intron_borders[0] and [1] refer to the next pair
                # del intron_borders[:2]

            # # if we encounter a STOP codon
            # if contig[n:n+3] in ['TAG','TAA','TGA']:
                # # create new exon and CDS features
                # new_exon = new_feature_from_template(g, 'exon')
                # new_CDS  = new_feature_from_template(g, 'CDS')

                # # # update new_exon properties
                # new_exon.start = exon_start
                # new_exon.end   = n+3
                # pattern        = "exon%02d" % exon_number 
                # new_exon.attributes['ID'][0] = new_exon.attributes['ID'][0].replace('exon01', pattern)
                # new_exon.id = new_exon.attributes['ID'][0]

                # # update new_CDS properties
                # new_CDS.start = exon_start
                # new_CDS.end   = n+3
                # new_CDS.frame = str(phase)
                # pattern       = "CDS%02d" % exon_number 
                # new_CDS.attributes['ID'][0] = new_CDS.attributes['ID'][0].replace('CDS01', pattern)
                # new_CDS.id = new_CDS.attributes['ID'][0]

                # # add new exon to new_features
                # new_features.extend([new_exon,new_CDS])
                
                # # edit gene start and end coordinates
                # new_gene.start = g.start                    # we assume start of gene feature is always correct
                # new_gene.end   = n+3
                # # add new gene to new_features
                # new_features.append(new_gene)

                # # we have reached the end of the gene,
                # # stop traversing the sequence
                # break

            # else:
                # n = n+3 
             
        # # delete the original + strand exon and CDS features
        # for f in db.children(g, featuretype=('exon','CDS')):
            # db.delete(f)
        # # delete the original + strand gene feature
        # db.delete(g)



# # update db with updated exon and CDS features
# db.update(new_features)


# # print updated gff record in logical order
# for f in db.all_features(order_by=('seqid','start','attributes'), limit='tig00000012:1-400000'):
#     if f.featuretype == 'gene':
#         print('\n', end='')
#         print(f)
#     else:
#         print(f)
