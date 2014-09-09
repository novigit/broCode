#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;
use Pod::Usage;

=head1 USAGE

    filterAssemblyLengthCov.pl -i|in [fasta-in] -l|length [length-cutoff] -c|cov [coverage-cutoff] -m|mode [AND|OR] > [fasta-out]

=head1 SYNOPSIS

    Takes contigs or scaffolds from an assembly, and selects those that meet the length cutoff and / or the cutoff, depending on the mode.

=head1 EXAMPLE

    filterAssemblyLengthCov.pl -i contigs.fasta -l 1000 -c 300 -m AND > contigs-filter.fasta
    Keeps all contigs longer than 1000 bp and coverage higher than 300x

    filterAssemblyLengthCov.pl -i contigs.fasta -l 1000 -c 300 -m OR > contigs-filter.fasta
    Keeps all contigs longer than 1000 bp or  coverage higher than 300x

=head1 COMMENTS

    Assumes that fasta header looks like ">NODE_[x]_length_[n]_cov_[n]"

=head1 TODO

    Modify it so that it can filter solely by length or coverage alone
    If filter by length alone, make it less dependent on format of the fasta header

=head1 AUTHOR

    Joran Martijn (joran.martijn@icm.uu.se)

=cut

pod2usage (-message => "not enough arguments") unless @ARGV >= 8;

my ($fa, $length_cutoff, $cov_cutoff, $mode);
GetOptions("i|in=s"     => \$fa,
	   "l|length=s" => \$length_cutoff,
	   "c|cov=s"    => \$cov_cutoff,
	   "m|mode=s"   => \$mode);

my $in  = Bio::SeqIO->new(-format => 'fasta',
			  -file   => $fa);
my $out = Bio::SeqIO->new(-format => 'fasta',
			  -fh     => \*STDOUT);

# Use ->next_seq to loop over the fasta file
while (my $seq = $in->next_seq) {

    # Use ->id() to retrieve fastaheader
    my $id = $seq->id;
    my @line = split "_", $id;
    my ($length, $cov) = ($line[3], $line[5]);

    # Use ->write_seq() to write sequences  ($seq) with sufficient length and cov to $out object
    if    ($mode eq "OR" ) {
	$out->write_seq($seq) if ($length > $length_cutoff || $cov > $cov_cutoff);
    }
    elsif ($mode eq "AND") {
	$out->write_seq($seq) if ($length > $length_cutoff && $cov > $cov_cutoff);
    }
}
