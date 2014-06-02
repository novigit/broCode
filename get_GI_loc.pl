#!/usr/bin/perl -w
# By Joran Martijn

use strict;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

=head1 USAGE

get_GI_loc.pl -i|--in [.gbk] -l|--giList [gi.tsv] -d|--distance [integer] > [.tsv]

=head1 COMMENTS

The gi's must be in a .tsv file looking something like this:

    flgB           347758455
    flgC           347758454
    flgD           347758481
[  PSIBLAST  ]       [GI]
[ QUERY GENE ]   

The outfile will be a .tsv file with a couple of different columns:

[GI]    [PSIBLAST QUERY GENE]    [LOCATION]     [GI GENE LENGTH]    [GI GBK ANNOTATION]

   
=head1 SYNOPSIS

Basically what this script does is parsing a genbank file while checking for GI numbers. If it identifies a GI number that you have provided, it will print out the GI, the gene which you used as query to find these GI's with a blast, the GI's location, the gene length, and the GI's gene annotation if available.

In addition, if the distance between two GI hits is greater than the threshold you provide with -d , it will print out '-------------'. This makes it easier to identify genes that are close to one another in the genome and potentially are in the same operon.

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se) with massive help from Lionel Guy (lionel.guy@icm.uu.se)

=cut

pod2usage (-message => "ERROR: NOT ENOUGH ARGUMENTS") unless @ARGV >= 6;

my ($gbk,$gi_list,$distance);
GetOptions ('i|in=s'       => \$gbk,
	    'l|giList=s'   => \$gi_list,
	    'd|distance=s' => \$distance);

my $in = Bio::SeqIO->new(-file   => $gbk,
			 -format => 'genbank');

open GI_LIST, '<', $gi_list;
my %hash;
while (<GI_LIST>) {
    chomp;
    my ($query_gene, $gi, $taxon) = split "\t", $_;
    # print $gene, "\t", $gi, "\n";
    $hash{$gi} = $query_gene;
}
close GI_LIST;

while (my $seq = $in->next_seq) {
    print "=" x 10, " ", $seq->id, " ", "=" x 10, "\n";
    my $last_end = 0;
    for my $ft ($seq->get_SeqFeatures){
	next unless ($ft->primary_tag eq 'CDS');
	next unless ($ft->has_tag('db_xref'));
	my @xrefs = $ft->get_tag_values('db_xref');
	foreach my $xref (@xrefs) {
	    next unless $xref =~ m/^GI:/;
	    # print $xref, "\n";
	    my ($gi) = $xref =~ m/GI:(\d+)/;
	    # print $gi, "\n";
	    if (exists $hash{$gi}) {
		print "--------------\n" if ($ft->start - $last_end > $distance);
		print $gi, "\t", $hash{$gi}, "\t", $ft->start(), "..", $ft->end(), "\t", $ft->length;
		print "\t", $ft->get_tag_values('gene') if ($ft->has_tag('gene'));
		print "\n";
		$last_end = $ft->end;
	    }
	}	
    }
}
