from Bio import SeqIO
from numpy import average
import argparse
import sys
import pysam
import textwrap

parser = argparse.ArgumentParser(
    description=textwrap.dedent(
        "get error rates per read from SAM file"
    ),
    epilog=textwrap.dedent("__doc__"),
    formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument("-i", "--input", required=True, help="SAM file")
parser.add_argument("-o", "--output", help="tab separated file with error rates", default=sys.stdout)
args = parser.parse_args()


cigar_tags = {'BAM_CMATCH':0,
                'BAM_CINS':1,
                'BAM_CDEL':2,
                'BAM_CREF_SKIP':3,
                'BAM_CSOFT_CLIP':4,
                'BAM_CHARD_CLIP':5,
                'BAM_CPAD':6,
                'BAM_CEQUAL':7,
                'BAM_CDIFF':8,
                'BAM_CBACK':9,
                'NM':10}


def main(infile, outfile):
    if isinstance(outfile, str):
        outfile = open(outfile, 'w')
    print("\t".join(["query", "alignment_length", "error_rate", "matches", "substitutions", "insertions", "deletions"]),
                        file=outfile)
    error_rates = []
    substitutions_overall = 0
    insertions_overall = 0
    deletions_overall = 0
    for ar in pysam.AlignmentFile(infile, 'r'):
        cigar_stats = ar.get_cigar_stats()[0]
        substitutions = cigar_stats[cigar_tags['NM']] - \
                (cigar_stats[cigar_tags['BAM_CINS']] + cigar_stats[cigar_tags['BAM_CDEL']])
        print("\t".join([ar.qname,
                        str(ar.query_alignment_length),
                        str(cigar_stats[cigar_tags['NM']] / ar.query_alignment_length),
                        str(cigar_stats[cigar_tags['BAM_CMATCH']] - substitutions),
                        str(substitutions),
                        str(cigar_stats[cigar_tags['BAM_CINS']]),
                        str(cigar_stats[cigar_tags['BAM_CDEL']])
                        ]), file=outfile)
        error_rates.append(cigar_stats[cigar_tags['NM']] / ar.query_alignment_length)
        substitutions_overall += substitutions
        insertions_overall += cigar_stats[cigar_tags['BAM_CINS']]
        deletions_overall += cigar_stats[cigar_tags['BAM_CDEL']]
    sum_errors = sum([substitutions_overall, insertions_overall, deletions_overall])
    print("average error rate in file %s was %0.4f%%" % (infile, average(error_rates)* 100) , file=sys.stderr)
    print("overall number of substitutions, insertions and deletions was %i (%0.4f%%), %i (%0.4f%%) and %i (%0.4f%%), respectively" %
                    (substitutions_overall, substitutions_overall / sum_errors * 100,
                    insertions_overall, insertions_overall / sum_errors * 100,
                     deletions_overall, deletions_overall / sum_errors * 100), file=sys.stderr)

if __name__ == '__main__':
    main(args.input, args.output)
