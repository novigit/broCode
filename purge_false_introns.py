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
# introns = db.features_of_type(featuretype='intron', limit='tig00000492:1-2167353')
introns = db.features_of_type(featuretype='intron')
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


# generate new feature objects from template
def new_feature_from_template(feat, ftype, strand):
    # determine attributes of template
    if ftype == 'gene':
        name = re.sub('ctg[0-9]{3}\.', '', feat.id) 
        attrs = 'ID='+feat.id+';Name='+name
    elif ftype =='mRNA':
        mrna_id = feat.id.replace('gene','mRNA')
        parent = feat.id
        name = re.sub('ctg[0-9]{3}\.', '', mrna_id) 
        attrs = 'ID='+mrna_id+';Parent='+parent+';Name='+name
    elif ftype =='intron':
        attrs = feat.attributes
    else:
        feat_id = feat.id.replace('gene','mRNA')
        attrs = 'ID='+feat_id+'.'+ftype+'01;Parent='+feat_id
    # define template feature
    feature_template = gffutils.feature.Feature(
            seqid=feat.seqid, 
            source='gffutils', 
            featuretype=ftype, 
            start=feat.start, 
            end=feat.end, 
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
    elif strand == '-':
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
    # sort selected orfs by start coordinate
    selected_orfs.sort(key=lambda x : x[0])
    return selected_orfs
    

# splice a sequence region with known intron coordinates
def splice_region(region, introns):
    '''
    Take a sequence string and a list of lists containing
    intron borders and splice out the introns from the 
    sequence string
    '''
    # -1 from intron start (='exon' end)
    # +1 from intron end   (='exon' start)
    introns = [ [b[0]-1 , b[1]+1] for b in introns ]
        
    # flatten list and convert to 0-indexing
    f = [ x-1 for b in introns for x in b ]
    # add start and end (0-indexed)
    f.insert(0, 0)
    f.append(len(region)-1)

    # repack list
    borders = [ [f[i],f[i+1]] for i in range(0, len(f), 2) ]
    # construct spliced region
    spliced_region = ''
    for b in borders:
        spliced_region += region[ b[0]:b[1]+1 ]
    return spliced_region


all_gene_ids = [ f.id for f in db.all_features(featuretype='gene') ]

new_features = []
# loop over genes with false introns
for g in genes_with_false_introns:

    # print()
    # print(g.id)

    # define region and strand on which we have to operate
    ## get start coordinate
    region_start = g.start
    ## get end coordinate (i.e. start of next gene - 1)
    gene_number = int( g.id.split('gene')[1] )
    gene_contig = str( g.id.split('gene')[0] )
    next_gene_id = gene_contig + "gene%04d" % (gene_number + 1)
    if next_gene_id in all_gene_ids:
        region_end = db[next_gene_id].start - 1
    else:
        region_end = g.end + 1
    ## extract region sequence
    ### note that slicing is 0-indexed, and gff coordinates 1-indexed
    contig = fasta[g.seqid].seq
    region_seq = str( contig[region_start-1:region_end] )
    ## get strand
    strand = g.strand

    # determine if there are any true introns in this gene
    ## get true intron borders
    introns = list( db.children(g, featuretype='intron', order_by=('seqid','start','attributes')) )
    true_introns = list( filter(lambda x : x.attributes['call'][0] == 'true', introns) )
    intron_borders = [ [t.start, t.end] for t in true_introns ]

    ## delete false introns
    false_introns = filter(lambda x : x.attributes['call'][0] == 'false', introns)
    db.delete(false_introns)

    # get gene ID
    gene_ID = g.attributes['ID'][0]

    # if there are no introns
    if len(true_introns) == 0:

        # apply naive ORF finder
        new_orfs = get_non_overlapping_orfs(region_seq, strand)
        n=1 # orf counter
        for orf in new_orfs:
            # create new template gene and mRNA feature
            new_gene = new_feature_from_template(g, 'gene', strand)
            new_mRNA = new_feature_from_template(g, 'mRNA', strand)
            new_exon = new_feature_from_template(g, 'exon', strand)
            new_CDS  = new_feature_from_template(g, 'CDS',  strand)
            # adjust coordinates
            new_gene.start = new_mRNA.start = new_exon.start = new_CDS.start = region_start + orf[0]
            new_gene.end   = new_mRNA.end   = new_exon.end   = new_CDS.end   = region_start + orf[1] - 1
            # update name only if multiple orfs
            if len(new_orfs) > 1:
                # update gene ID and Name
                new_gene_ID = re.sub('([0-9]{4})', r'\1-' + str(n), gene_ID)
                new_gene.id = new_gene.attributes['ID'][0] = new_gene_ID
                new_gene.attributes['Name'][0] = re.sub('ctg[0-9]{3}\.', '', new_gene_ID)
                # update mRNA name
                new_mRNA_ID = re.sub('gene', 'mRNA', new_gene_ID)
                new_mRNA.id = new_mRNA.attributes['ID'][0] = new_mRNA_ID
                new_mRNA.attributes['Parent'][0] = new_gene_ID
                new_mRNA.attributes['Name'][0] = re.sub('ctg[0-9]{3}\.', '', new_mRNA_ID)
                # update exon name
                new_exon_ID = new_mRNA_ID + '.exon01'
                new_exon.id = new_exon.attributes['ID'][0] = new_exon_ID
                new_exon.attributes['Parent'][0] = new_mRNA_ID
                # update CDS name
                new_CDS_ID = new_mRNA_ID + '.CDS01'
                new_CDS.id = new_CDS.attributes['ID'][0] = new_CDS_ID
                new_CDS.attributes['Parent'][0] = new_mRNA_ID

            # print(new_gene)
            # print(new_mRNA)
            # print(new_exon)
            # print(new_CDS)

            # update orf counter
            n = n+1

            # add new features to list
            new_features.extend( [new_gene,new_mRNA,new_exon,new_CDS] )

    # if there are introns, however..
    elif len(true_introns) != 0:

        # adjust intron coordinates so they are
        # relative to the start of the sequence region
        rel_intron_borders = [ [start-region_start+1,stop-region_start+1] for start, stop in intron_borders ]

        # splice out introns from sequence region
        spliced_region_seq = splice_region(region_seq, rel_intron_borders)

        # calculate intron junction coordinates in spliced sequence region
        # and intron length
        junctions = []
        accum_intron_len = 0
        for i in rel_intron_borders:
            junction = i[0] - accum_intron_len
            intron_len = i[1] - i[0] + 1
            junctions.append([ junction,intron_len ])
            accum_intron_len += intron_len

        # apply naive ORF finder
        new_orfs = get_non_overlapping_orfs(spliced_region_seq, strand)
        # convert to list of lists so I can update orf coordinates
        # (tuples are immutable)
        upd_orfs = [ list(orf) for orf in new_orfs ]

        # re-insert introns
        # (i.e. update orf start and end coordinates)
        accum_intron_len = 0
        # compare all-vs-all orfs vs intron-junctions
        for o in upd_orfs:
            # shift orf positions by accumulated intron lengths
            o[0] += accum_intron_len
            o[1] += accum_intron_len
            for j in junctions:
                i,intron_len = j
                # if intron junction is situated before the orf
                if i < o[0] < o[1]:
                    o[0] += intron_len
                    o[1] += intron_len
                    junctions.remove(j) # intron is 'used', no longer consider it
                    accum_intron_len += intron_len
                # if intron junction is situated within the orf
                if o[0] < i < o[1]:
                    o[1] += intron_len
                    junctions.remove(j) # intron is 'used', no longer consider it
                    accum_intron_len += intron_len    
                # if intron junction is situated after  the orf    
                if o[0] < o[1] < i:
                    pass

        # convert updated ORF coordinates back to genome space
        upd_orfs  = [ [start+region_start,stop+region_start-1,strand,desc] for start,stop,strand,desc in upd_orfs ]

        # transfer updated ORF coordinates and exon coordinates to new features
        # all-vs-all orf vs intron comparison
        n=1 # orf counter
        for o in upd_orfs:
            new_gene = new_feature_from_template(g, 'gene', strand)
            new_mRNA = new_feature_from_template(g, 'mRNA', strand)
            new_gene.start = new_mRNA.start = exon_start = o[0]
            new_gene.end   = new_mRNA.end   = exon_end   = o[1]

            # gene_ID = 'ctg012.gene0001'
            new_mRNA_ID = re.sub('gene','mRNA', gene_ID)
            # print(new_mRNA_ID)

            # update gene names only if there are multiple new orfs
            if len(upd_orfs) > 1:
                # update gene ID and Name
                new_gene_ID = re.sub('([0-9]{4})', r'\1-' + str(n), gene_ID)
                new_gene.id = new_gene.attributes['ID'][0] = new_gene_ID
                new_gene.attributes['Name'][0] = re.sub('ctg[0-9]{3}\.', '', new_gene_ID)
                # update mRNA ID and Name and Parent
                new_mRNA_ID = re.sub('gene', 'mRNA', new_gene_ID)
                new_mRNA.id = new_mRNA.attributes['ID'][0] = new_mRNA_ID
                new_mRNA.attributes['Name'][0] = re.sub('ctg[0-9]{3}\.', '', new_mRNA_ID)
                new_mRNA.attributes['Parent'][0] = new_gene_ID


            # print(new_gene)
            # print(new_mRNA)

            # m=1 # exon counter
            # for i in intron_borders:
            #     # fish out intron_object
            #     intron = true_introns[m-1]
            #     # print(intron)
            m=1 # exon counter
            for intron in true_introns:
                i = [ intron.start, intron.end ]
                # print(i[0])
                if i[0] < o[0] < o[1]:
                    # intron_borders.remove(i)
                    true_introns.remove(intron)
                    db.delete(intron)
                if o[0] < i[0] < o[1]:
                    # create new exon and CDS
                    new_exon       = new_feature_from_template(g, 'exon', strand)
                    new_CDS        = new_feature_from_template(g, 'CDS', strand)
                    ## update exon coordinates
                    new_exon.start = new_CDS.start = exon_start
                    new_exon.end   = new_CDS.end   = i[0] - 1

                    ## update exon ID and Parent
                    new_exon_ID    = new_mRNA_ID + ".exon%02d" % m
                    new_exon.id    = new_exon.attributes['ID'][0] = new_exon_ID
                    new_exon.attributes['Parent'][0] = new_mRNA_ID

                    ## update CDS ID and Parent
                    new_CDS_ID = new_mRNA_ID +  ".CDS%02d"  % m
                    new_CDS.id = new_CDS.attributes['ID'][0] = new_CDS_ID
                    new_CDS.attributes['Parent'][0] = new_mRNA_ID

                    # create new corresponding intron

                    # print()
                    # print(intron)
                    # print()

                    new_intron    = new_feature_from_template(intron, 'intron', strand)
                    new_intron_ID = new_mRNA_ID +  ".intron%02d"  % m
                    new_intron.id = new_intron.attributes['ID'][0] = new_intron_ID
                    new_intron.attributes['Parent'][0] = new_mRNA_ID

                    # print(new_exon)
                    # print(new_CDS)
                    # print(new_intron)

                    new_features.extend([ new_exon,new_CDS,new_intron ])
                    
                    # move exon start downstream
                    exon_start = i[1] + 1
                    
                    # update exon counter
                    m = m+1

                    # intron is 'used', no longer consider it
                    # intron_borders.remove(i)
                    true_introns.remove(intron)
                    db.delete(intron)

                if o[0] < o[1] < i[0]:
                   pass

            # update exon & CDS features
            ## if orf has no introns
            ## or if last exon of the orf
            new_exon = new_feature_from_template(g, 'exon', strand)
            new_CDS  = new_feature_from_template(g, 'CDS', strand)
            new_exon.start = new_CDS.start = exon_start
            new_exon.end   = new_CDS.end   = exon_end
            # update exon ID and Parent
            new_exon_ID = new_mRNA_ID + ".exon%02d" % m
            new_exon.id = new_exon.attributes['ID'][0] = new_exon_ID
            new_exon.attributes['Parent'][0] = new_mRNA_ID
            # update CDS ID and Parent
            new_CDS_ID = new_mRNA_ID +  ".CDS%02d"  % m
            new_CDS.id = new_CDS.attributes['ID'][0] = new_CDS_ID
            new_CDS.attributes['Parent'][0] = new_mRNA_ID

            # print(new_exon)
            # print(new_CDS)

            new_features.extend([ new_gene, new_mRNA, new_exon, new_CDS ])

            # update orf counter
            n = n+1


    # delete the original + strand exon and CDS features
    # print('TO BE DELETED')
    for f in db.children(g, featuretype=('exon','intron','CDS','mRNA')):
        # if f.featuretype == 'intron':
        #     print(f)
        db.delete(f)
    # delete the original + strand gene feature
    db.delete(g)

# print(len(new_features))
# print(new_features)

# for f in new_features:
#     print(f)

# update db with updated exon and CDS features
db.update(new_features)


# print updated gff record in logical order
for f in db.all_features(order_by=('seqid','start','attributes')):#, limit='tig00000012:1-400000'):
    if f.featuretype == 'gene':
        print()
        print(f)
    else:
        print(f)

