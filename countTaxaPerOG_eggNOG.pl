#!/usr/bin/perl -w
use strict;

=head1 NAME

countTaxaPerOG_eggNOG.pl - Create a count-table from an eggNOG <eggNOG>.member.tsv

=head1 USAGE

countTaxaPerOG_eggNOG.pl -i <eggNOG>.member.tsv > <count-table>.tsv

=head1 DESCRIPTION

Parses a <eggNOG>.member.tsv to count the number of times each taxon occurs in each OG in this eggNOG. Output in a tab separated value table

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)

=cut

# libraries
use Getopt::Long;
use Pod::Usage;

# state options
my ($input, $help);
GetOptions(
    'i|input=s' => \$input,
);
pod2usage (-msg => "Not enough arguments") unless $input;


# datastructures
my %count;
my %taxa;

# parse file
open FILE, '<', $input;
while (my $line = <FILE>) {

    # extract OG id and taxids
    my ($og,$content) = ( split '\t', $line )[1,5];
    my @seqids = split ',', $content;
    foreach my $seqid (@seqids) {
	my ($taxid) = $seqid =~ m/(\d+)\..+/;
	# store and update count 
	$count{$og}{$taxid}++;
	$taxa{$taxid}++;
    }
}
close FILE;

# print table header (taxids)
my @taxa = sort keys %taxa;
print "eggNOG", "\t", join("\t", @taxa), "\n";

# print rows (each row is one OG, each cell is a count)
foreach my $og (sort keys %count) {
    print $og, "\t";
    foreach my $taxid (@taxa) {
    	if ( $count{$og}{$taxid}) {print $count{$og}{$taxid}}
    	else                      {print "0"}
    	print "\t";
    }
    print "\n";
}
