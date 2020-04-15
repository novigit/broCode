#!/usr/bin/env python
import re
import pdb


def make_gf2gname_and_gflist(gf_file):
    gf2name = {}
    gflist = []
    with open(gf_file) as f:
        for line in f:
            gf, gname = line.rstrip('\n').split('\t')
            gf2name[gf] = gname
            gflist.append(gf)
    return gf2name, gflist


# map all non-fusion gfs to their nog
def make_gf2nog2arnogs_simple_gfs(gf2name, emapper_file):
    gf2nog2arnogs = {}
    gfs_of_interest = gf2name.keys()
    with open(emapper_file) as f:
        for line in f:
            # skip line if it does not contain a gf of interest
            # explicitly stating continue saves an indent level
            if not any(gf in line for gf in gfs_of_interest):
                continue
            # get gf, nog and arnogs
            gf = re.search('(.{5})@eurNOG', line).group(1)
            nog = re.search('(\\w+)@NOG', line).group(1)
            arnogs = set(re.findall('(\\w+)@arNOG', line))
            # if gf2nog2arnogs[gf] does not exist declare it
            if gf not in gf2nog2arnogs.keys():
                gf2nog2arnogs[gf] = {}
                gf2nog2arnogs[gf][nog] = arnogs
            # else if it does exist update the arnogs
            else:
                # gf2nog2arnogs[gf][nog] returns the
                # set of arnogs encountered in
                # previous lines with this gf
                gf2nog2arnogs[gf][nog].update(arnogs)
    return gf2nog2arnogs


def make_gf2nog2arnogs_fusion_gfs(gf2name, arcog_file):
    # get gf2nog mapping
    # from gf2name dict
    gf2cog = {}
    for gf in gf2name.keys():
        if '_COG' in gf:
            # in arCOGdef.tab file COG identifiers have an extra 0
            cog = re.sub('\\w{5}_COG', 'COG0', gf)
            gf2cog[gf] = cog

    # get nog2arnogs mapping
    # from arcog file
    cog2arcogs = {}
    with open(arcog_file) as f:
        for line in f:
            # skip lines that do not have COGs of interest
            if not any(cog in line for cog in gf2cog.values()):
                continue
            arcog, cog = [line.rstrip('\n').split('\t')[i] for i in (0, 4)]
            arcogs = cog2arcogs.get(cog, [])
            arcogs.append(arcog)
            cog2arcogs[cog] = arcogs

    # compile gf2nog2arnogs from
    # gf2nog and cog2arcogs
    gf2nog2arnogs = {}
    for gf, nog in gf2cog.items():
        # get list of arnogs associated with nog
        # if nog doesnt exist in cogs2arcogs
        # (the nog has no associated arcogs)
        # return empty ''
        arnogs = cog2arcogs.get(nog, '')
        gf2nog2arnogs[gf] = {}
        gf2nog2arnogs[gf][nog.replace('COG0', 'COG')] = set(arnogs)

    return gf2nog2arnogs


def merge_two_dicts(dict1, dict2):
    dict3 = dict1.copy()
    dict3.update(dict2)
    return dict3


def get_eurnog2funcat2annot(gf2name, eurnogs_file):
    # get eurnog2funcat2annot mapping
    eurnog2funcat2annot = {}
    # trim '_COGxxxx' because that notation doesnt exist in the eurnogs_file
    eurnogs_of_interest = [re.sub('_COG[0-9]{4}', '', gf) for gf in gf2name.keys()]
    with open(eurnogs_file) as f:
        for line in f:
            # skip line if it does not contain a gf of interest
            # explicitly stating continue saves an indent level
            if not any(eurnog in line for eurnog in eurnogs_of_interest):
                continue
            fields = line.replace('ENOG41', '').rstrip('\n').split('\t')
            eurnog, funcat, annot = [fields[i] for i in (1, 4, 5)]
            # extract eurnog corresponding gf from gf2name.keys()
            eurnog2funcat2annot[eurnog] = (funcat, annot)
    return eurnog2funcat2annot


def get_nog2funcat2annot(gf2nog2arnogs_all, nogs_file, cogs_file):
    # compile nogs of interest 
    nogs_of_interest = set([list(d.keys())[0] for d in gf2nog2arnogs_all.values()])

    # declare nog2funcat2annot mapping
    nog2funcat2annot = {}

    # get nog annotations from the eggNOG annotation file
    # for those NOGs that are not COGs
    nog_nogs = [nog for nog in nogs_of_interest if 'COG' not in nog]
    with open(nogs_file) as nf:
        for line in nf:
            # skip line if it does not contain a gf of interest
            # explicitly stating continue saves an indent level
            if not any(nog in line for nog in nog_nogs):
                continue
            fields = line.replace('ENOG41', '').rstrip('\n').split('\t')
            nog, funcat, annot = [fields[i] for i in (1, 4, 5)]
            # extract nog corresponding gf from gf2name.keys()
            nog2funcat2annot[nog] = (funcat, annot)

    # get nog annotations from the COG annotation file
    # for those NOGs that are COGs
    cog_nogs = [nog for nog in nogs_of_interest if 'COG' in nog]
    with open(cogs_file, encoding="ISO-8859-1") as cf:
        for line in cf:
            # skip line if it does not contain a gf of interest
            # explicitly stating continue saves an indent level
            if not any(cog in line for cog in cog_nogs):
                continue
            cog, funcat, annot = line.replace('ENOG41', '').rstrip('\n').split('\t')
            # extract nog corresponding gf from gf2name.keys()
            nog2funcat2annot[cog] = (funcat, annot)

    return nog2funcat2annot


def get_arnogs2funcat2annot(gf2nog2arnogs_all, arnogs_file, arcog_file):
    # compile set of arnogs of interest
    arnogs_of_interest = set([item for subdict in gf2nog2arnogs_all.values() for valueset in subdict.values() for item in valueset])

    # declare arnog2annot mapping
    arnog2annot = {}

    # get arNOG (those without arCOG equivalents) from the eggnog annotation file
    with open(arnogs_file) as f:
        for line in f:
            # skip line if it does not contain a gf of interest
            if not any(arnog in line for arnog in arnogs_of_interest if 'arCOG' not in arnog):
                continue
            arnog, annot = [line.replace('ENOG41', '').rstrip('\n').split('\t')[i] for i in (1, 5)]
            arnog2annot[arnog] = annot

    # get arCOG annotations from the arCOG annotation file
    with open(arcog_file) as f:
        for line in f:
            # skip line if it does not contain a gf of interest
            if not any(arcog in line for arcog in arnogs_of_interest if 'arCOG' in arcog):
                continue
            arcog, annot = [line.rstrip('\n').split('\t')[i] for i in (0, 3)]
            arnog2annot[arcog] = annot

    return arnog2annot


def print_table_sorted_by_gflist(gflist, gf2name, gf2nog2arnogs, eurnog2funcat2annot, nog2funcat2annot, arnog2annot):
    for gf in gflist:

        # gene name, nog and arnog ids
        gname = gf2name.get(gf)

        # eurnog id, funcat and annot
        eurnog = re.sub('_COG[0-9]{4}', '', gf)
        (eurnog_funcat, eurnog_annot) = eurnog2funcat2annot.get(eurnog)

        # nog id, funcat and annot
        nog = list(gf2nog2arnogs.get(gf).keys())[0]
        (nog_funcat, nog_annot) = nog2funcat2annot.get(nog)

        # arnog ids and their annots
        arnogs = gf2nog2arnogs.get(gf).get(nog)
        arnogs_annots = [arnog2annot.get(arnog, 'NO_ARNOG_ANNOTATION_FOUND') for arnog in arnogs]

        # print everything
        try:
            print(gf + '\t' + gname + '\t' + eurnog_funcat + '\t' + eurnog_annot + '\t' + nog + '\t' + nog_funcat + '\t' + nog_annot + '\t' + ' | '.join(arnogs) + '\t' + ' | '.join(arnogs_annots))
        except:
            pdb.post_mortem()

# get gene families of interest from file
gf2name, gflist = make_gf2gname_and_gflist('gfs_of_interest.csv')

# get gene family-to-nog-to-arnog mapping for simple non-fusion gfs
gf2nog2arnogs_simple_gfs = make_gf2nog2arnogs_simple_gfs(gf2name, 'methanotecta.emapper.annotations.all')

# get gene family-to-nog-to-arnog mapping for fusion associated gfs
gf2nog2arnogs_fusion_gfs = make_gf2nog2arnogs_fusion_gfs(gf2name, 'ar14.arCOGdef.tab')

# compile overarching gf2nog2arnogs_all dictionary
gf2nog2arnogs_all = merge_two_dicts(gf2nog2arnogs_simple_gfs, gf2nog2arnogs_fusion_gfs)

# get gene family-to-funcat-to-annotation mapping for eurnogs
eurnog2funcat2annot = get_eurnog2funcat2annot(gf2name, 'eurNOG.annotations.tsv')

# get gene family-to-funcat-to-annotation mapping for nogs
nog2funcat2annot = get_nog2funcat2annot(gf2nog2arnogs_all, 'NOG.annotations.tsv', 'cognames2003-2014.tab')

# get arnog-to-annotation mapping
arnog2annot = get_arnogs2funcat2annot(gf2nog2arnogs_all, 'arNOG.annotations.tsv', 'ar14.arCOGdef.tab')

# print entire table
print_table_sorted_by_gflist(gflist, gf2name, gf2nog2arnogs_all, eurnog2funcat2annot, nog2funcat2annot, arnog2annot)

# pdb.set_trace()
