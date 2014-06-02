#!/usr/bin/perl -w
use strict;

=head1 NAME 

assign-genes2domains.pl

=head1 USAGE

assign-genes2domains.pl [in]megan.paths [out]domain_count.tsv

in : MEGAN .paths  file that lists for each blast-query its taxonomic path with confidence values

out: Tab seperated file that lists for each contig how many blast-queries (genes) belong to which taxonomic domain or phylum 

=head1 COMMENTS

-This script works only with MEGAN path files in which the blast-queries have the following format:
'NODE_x_length_x_cov_x_[contig-count]_[gene-count]'

-The reported number of Bacterial and Archaeal genes does not include the amount of Alphaproteobacterial and Methanobacterial genes, respectively

-The output of this script can be used as input in the assign-contigs2domains.pl script


=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)
   
=cut

die "usage:assign-genes2domains.pl [megan.paths] [out.tsv]", "\n" unless @ARGV == 2;

my ($in, $out) = @ARGV;

my $id;
my %count;
my @taxa = ('Archaea',
	    'Methanobacteria',
	    'Bacteria',
	    'Alphaproteobacteria',
	    'Eukaryota', 
	    'Viruses', 
	    'No hits');

open IN, "<$in";
while (<IN>) {
    chomp;
    my ($node) = split ';', $_;
    my @node = split '_', $node;
    my $id = join '_', @node[0..5];
    print $id, "\n";

    foreach my $taxon (@taxa) {
	$count{$id}{$taxon}++ if (/$taxon/);
    }

    $count{$id}{'Archaea'}--  if (/Methanobacteria/);
    $count{$id}{'Bacteria'}-- if (/Alphaproteobacteria/);
}
close IN;

open OUT, ">$out";
print OUT "NODE", "\t";
foreach my $taxon (@taxa) {print OUT $taxon, "\t"}
print OUT "\n"; 

foreach my $node (keys %count ) {
    print OUT $node, "\t"; 
    foreach my $taxon (@taxa) {
	if (exists $count{$node}{$taxon})   {print OUT $count{$node}{$taxon}}
	else                                {print OUT "0"}
	print OUT "\t";
    }
    print OUT "\n";
}
close OUT;


