#!/usr/bin/perl -w

=head1 NAME

removeSequences.pl

=head1 SYNOPSIS

This script will remove contigs from a FASTA file that you specify. 

=head1 COMMENTS

Used to be called 'select_contigs.pl' and 'removeContigs.pl'.
Currently only works with one query only. Needs to be updated.

Works with both nucleotide and protein fasta

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)

=cut

use strict;
use Bio::SeqIO;

die "usage: removeSequences.pl [fasta-in] [seq1] [seq2] [etc] > [fasta-out]\n" unless @ARGV >= 2;

my ($fa, @queries) = @ARGV;

my $in  = Bio::SeqIO->new(-format => 'fasta',
			  -file   => $fa);
my $out = Bio::SeqIO->new(-format => 'fasta',
			  -fh     => \*STDOUT);


while (my $seq = $in->next_seq) {
    my $id = $seq->display_id;
    # print $id, "\n";
    my $print = 1;
    foreach my $query (@queries) {
	$print = 0 if ($id =~ m/$query/);
    }
    $out->write_seq($seq) if $print;
}    
