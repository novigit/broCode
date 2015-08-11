#!/usr/bin/perl -w

=head1

Usage: filterFastaFromLength.pl -c -s min_size multi_fasta1 multi_fasta2 ...

min_size       Minimal size required for a contig to be included
multi_fastax   (Multi-)fasta file(s)
-c             Check unique contig names
-f             Format of the file (fasta, qual,...)


Function: filters several fasta or multi-fasta files and merge them in one
          single file. Additionnaly, can filter from contig size and check
          that contigs don't have the same id.

=cut

use strict;
use Bio::SeqIO;
use Getopt::Std;

our $opt_s;
our $opt_c = 1; # perform name check?
our $opt_f = 'fasta';
&getopts('cs:f:');

my %seen;
foreach (@ARGV){
    my $seq_in = Bio::SeqIO->new(
	-format => $opt_f,
	-file   => $_
    );
    my $seq_out = Bio::SeqIO->new(
	-format => $opt_f,
	-fh     => \*STDOUT
    );
    while (my $seq = $seq_in->next_seq) {
	my $id = $seq->id;
	# if check enabled, check that contig name is unique
	die ("Non-unique fasta record names\n") if ($opt_c && $seen{$id}++);
	if ( $seq->length >= $opt_s ){
	    $seq_out->write_seq($seq);
	}
    }
}
