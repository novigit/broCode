
from ETE3_Utils import *
from ETE3_styles import *
from misc_utils import *
from operator import gt, lt
import argparse
import os
import warnings

def get_counts(taxon_prefix, prefix_taxon):
    """
    Counts how many times a value appears in dictionary (as value) and
    returns another dictionary with the counts.
    :param taxon_prefix:
    :param prefix_taxon:
    :return:
    """
    taxon_copies = {}
    for taxon, prefix in taxon_prefix.items():
        taxon_copies[taxon] = len(prefix_taxon[prefix])
    return taxon_copies


def parse_args():
    ## Parse arguments
    parser = argparse.ArgumentParser(description="This script coverts a tree into pdf. It can root the tree, colour it and highlight repeated taxa."
                                                 "Right now the it only works with the option --prefixes. This option assumes that"
                                                 "the mapping file has two columns <prefix> and <clade> and that the outgroup file "
                                                 "contains prefixes. Prefixes will be assumed"
                                                 "to be the name of the organism when --duplicated is set. YOU NEED TO HAVE python_lib IN"
                                                 "YOUR PYTHONPATH.")

    parser.add_argument("-t", "--tree", required=True, help="Newick tree file")
    parser.add_argument("-p", "--pdf", required=True, help="Pdf name")
    parser.add_argument("-m", "--mapping", required=False, default=None, help="Mapping file. Tab separated")
    parser.add_argument("-o", "--outgroup", required=False, default=None, help="File including taxa to be used as outgroup. One name per line. "
                                                                               "If --duplicated is set it will exclude those taxa in the outgroup that are duplicated")
    parser.add_argument("--prefixes", required=False, default=False, action='store_true', help="Mapping file of prefixes. "
                                                                                               "It requires -m.")
    parser.add_argument("-d", "--duplicated", required=False, default=False, action='store_true', help="Highlight duplicated taxa."
                                                                                                       "It requires -m")
    parser.add_argument("-c", "--color", required=False, default="random", help="Color palette to use. By default random. \n"
                                                                                "Available palettes: \n\t%s" % ("\n\t".join(list(AVAIL_PALETTES))))

    args = parser.parse_args()

    return args


def main(tree_file, pdf, mapping_file, outgroup_file, prefix, duplicated, palette):

    # Read tree
    #tree = PhyloTree(tree_file)
    #tree = parse_newick(tree_file)
    try:
        tree = parse_newick(tree_file)
    except:
        warnings.warn("Opppss.. there is some problem with your newick file. I will try to read it but the support values might go bananas")
        tree = Tree(tree_file, format=1)

    set_node_style(tree, node_style_basic)

    if mapping_file is not None:

        if prefix:
            # Read mapping file of prefixes [prefix(taxa)   clade]
            taxon_clade, clade_taxon, taxon_prefix, prefix_taxon = read_prefix_map(tree, mapping_file)

        else:
            # Read mapping file of prefixes [leaves_name   clade  (optional:taxa)]
            print "Read from names-mapping file. TODO. Not implemented yet"

        #print len(taxon_clade), len(clade_taxon), len(taxon_prefix), len(prefix_taxon)

        # Set palette
        if palette in AVAIL_PALETTES:
            palette = AVAIL_PALETTES[palette]
        else:
            palette = assign_random_colors(list(clade_taxon.keys()))

        # Add features to tree
        initiate_feature(tree, taxon_clade, "clade")
        initiate_feature(tree, combine_dictionaries(taxon_clade, palette), "color")


    # Highlight marker genes and exclude them from outgroup if they are duplicated
    if duplicated:
        if prefix:
            # Add counts number
            initiate_feature(tree, get_counts(taxon_prefix, prefix_taxon), "counts")
            root_cond = (lt, "counts", 2)
            set_node_style(tree, node_style_highlight, leaves=True, condition=(gt, "counts", 1))
            #initiate_feature(tree, taxon_prefix, "organism")
        else:
            # TODO
            pass
    else:
        # TODO
        root_cond = None

    # Root tree
    if outgroup_file is not None:
        if prefix:
            # Read outgroup
            outgroup = []
            for line in open(outgroup_file):
                outgroup_prefix = line.strip()
                outgroup += prefix_taxon[outgroup_prefix]
            root_tree(tree, outgroup, root_cond)

    # Set styles and make pdf
    tree.ladderize(direction=1)
    title = os.path.splitext(os.path.basename(pdf))[0]
    ts = tree_style_basic(layout_node_color, title)
    #tree.show(tree_style=ts)
    tree.render(pdf, w=1500, units="px", tree_style=ts)

    # Write nexus file
    write_nexus(tree, os.path.splitext(pdf)[0] + ".nexus")


if __name__ == "__main__":
    args = parse_args()
    main(args.tree, args.pdf, args.mapping, args.outgroup,
         args.prefixes, args.duplicated, args.color)
