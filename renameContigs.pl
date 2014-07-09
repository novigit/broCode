#!/usr/bin/perl -w

=head1 USAGE

renameContigs.pl -f|fasta [.fasta] -t|table [.tsv] -d|digit [integer] > [out.fasta]

=head1 SYNOPSIS

Renames contigs from an assembly into NCBI submittable headers and easy prokka processing. 

    For example: '>NODE_52_length_456_cov_103.6' to '>contig001'

You can set the amount of digits to fit the number of contigs. 
For example if you have 10000 contigs you would like to have 5 digits.
But if you have 340 contigs you would only need 3 digits.

Also creates a tsv file or dictionary where you can see which contig name belongs to which assembly header.

For example:

    NODE_52_length_456_cov_103.6     contig001
    NODE_58_length_50231_cov_500.2   contig002
    NODE_68_length_236_cov_10000.1   contig003

=cut

use strict;
use Getopt::Long;
use Pod::Usage;

pod2usage (-message => "Not enough arguments") unless @ARGV >= 6;

my ($fasta, $table, $digit);
GetOptions('f|fasta=s' => \$fasta,
	   't|table=s'  => \$table,
	   'd|digit=s'  => \$digit);

my $count = 0;
my $format = '%0' . $digit . 'd';
my %translate;

open FASTA, "<$fasta";
open TABLE, ">$table";
while (<FASTA>){
    if (/^>/) {

	# write out fasta
	$count++;
	my $contig = 'contig' . sprintf("$format", $count);
	my ($header) = $_ =~ m/^>(.+)/;
	$_ =~ s/>$header/>$contig/;
	print $_;

	# write out table (dictionary)
	print TABLE $header, "\t", $contig, "\n";
    }
    else {print $_};
}
close TABLE;
close FASTA;
