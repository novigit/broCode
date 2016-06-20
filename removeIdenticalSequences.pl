#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;

# Based on id

my %unique;

my ($file) = @ARGV;
my $in     = Bio::SeqIO->new(-file => $file,    -format => "fasta");
my $out    = Bio::SeqIO->new(-fh   => \*STDOUT, -format => "fasta");

# loop over sequences
while(my $seq_obj = $in->next_seq) {

    # get header id
    my $id  = $seq_obj->display_id;
   
    # print sequence only if header does not exist in hash yet
    unless(exists($unique{$id})) {
	$out->write_seq($seq_obj);
	$unique{$id}++;
    }
}
