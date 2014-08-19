#!/usr/bin/perl

# Script will extract gene based on a identifier in the header, can be gi, organism name or whatever
# By Anders Lind, edited by Joran Martijn
# Opposite of removeSequences.pl

=head1 USAGE

selectSequences.pl -i [fasta-in] -l [list] -q [query(ies)] > [fasta-out]

Either use -l or -q

-q Query can be partial of sequence header

-l Query has to be exactly the same as sequence header

=cut 

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

pod2usage (-message => "Not enough arguments") unless @ARGV >= 4;

my ($fasta, $list, $queries);
GetOptions("i|fasta=s"   => \$fasta,
	   "l|list=s"    => \$list,
	   "q|queries=s" => \$queries);

my $in  = Bio::SeqIO->new(-format => 'fasta',
			  -file   => $fasta);
my $out = Bio::SeqIO->new(-format => 'fasta',
			  -fh     => \*STDOUT);

if ($queries) {
    
    my @queries = split(/,/, $queries);
    my $count = 0;

    while (my $seq = $in->next_seq) {
    	my $id = $seq->id;
    	$id .= $seq->desc;
	# print $id, "\n";
	
	foreach my $query (@queries){
	    # print $query, "\n";
	    if ($id =~ m/$query/){
	    	$count++;
	    	$out->write_seq($seq);
	    }
	}
    }

# Rewrite this part
    foreach my $query (@queries){
	print STDERR "Could not find any headers containing $query\n" if $count == 0;
	print STDERR "Found $count header(s) containing $query\n"     if $count > 0;
    }

}

if ($list) {
    my %hash;
    
    open LIST, "<$list";
    while (<LIST>) {
	chomp;
	my $query = $_;
	# print $query, "\n";
	$hash{$query}++;
    }
    close LIST;
    
    while (my $seq = $in->next_seq) {
    	my $id = $seq->id;
    	$id .= $seq->desc;
    	# print $id, "\n";
    	# $out->write_seq($seq) if (exists $hash{$id});
	
    	my $print = 0;
    	foreach my $query (keys %hash) {
    	    $print = 1 if ($id =~ m/$query/);
    	}
    	$out->write_seq($seq) if $print;
    }
}
