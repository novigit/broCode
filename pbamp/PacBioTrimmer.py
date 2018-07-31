#! /usr/bin/env python3
"""Trim PacBio reads with stretches of low quality"""
import argparse
from Bio import SeqIO
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", required=True,
                    help="input fastq file")
parser.add_argument("-t", "--threshold", required=False, default=18, type=int,
                    help="threshold for low quality stretches, the program looks \
                    for stretches that have a mean below this value and removes \
                    the reads, default 18")
parser.add_argument("-w", "--window", required=False, default=30, type=int,
                    help="sliding window size, default 30")
parser.add_argument("--print", dest='print', default=False, action='store_true',
                    help="print a png heatmap for the discarded reads, \
                    very(!) slow but might be interesting, especially for debugging")



args = parser.parse_args()


def containsLowQualityStretches(readname, quality, threshold, window):
    if any([score < threshold for score in quality.rolling(window).mean().dropna().tolist()]):
        print("discard %s" % readname)
        return True
    else:
        return False


def main(args):
    discarded_count = 0
    kept_count = 0
    discarded_qualities = []
    with open("%s_discarded.fastq" % args.input.replace('.fastq', ''), 'w') as discarded, \
        open("%s_kept.fastq" % args.input.replace('.fastq', ''), 'w') as kept:
        for s in SeqIO.parse(args.input, 'fastq'):
            quality = pd.Series(s.letter_annotations['phred_quality'])
            if containsLowQualityStretches(s.id, quality, args.threshold, args.window):
                SeqIO.write(s, discarded, "fastq")
                discarded_count += 1
                discarded_qualities.append(quality)
            else:
                SeqIO.write(s, kept, "fastq")
                kept_count += 1
    if args.print:
        df = pd.concat([pd.DataFrame(q) for q in discarded_qualities], join='outer')
        sns.heatmap(df)
        plt.savefig("discarded_qualities.png")
    print("##############################################")
    print()
    print("discarded %s reads of a total of %s" % (discarded_count,
                                                discarded_count + kept_count))

if __name__ == '__main__':
    main(args)
