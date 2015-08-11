#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;

# Based on id

my %hash;

my ($file) = @ARGV;
my $in     = Bio::SeqIO->new(-file => $file,    -format => "fasta");
my $out    = Bio::SeqIO->new(-fh   => \*STDOUT, -format => "fasta");

# loop over sequences
while(my $seq = $in->next_seq) {

    # get header id and length of seq
    my $id     = $seq->display_id;
    my $length = $seq->length;
    # print $id, "\t", $length, "\n";

    # if new id, store $seq and $length in hash
    unless (exists ($hash{$id})) {
	# $hash{$id}++;
    	$hash{$id}{'seq'}    = $seq;
    	$hash{$id}{'length'} = $length;
    }

    # if identical id, overwrite $seq only if sequence is larger
    elsif ($length > $hash{$id}{'length'}) {
	# print "i get here", "\n";
    	$hash{$id}{'seq'} = $seq;
    }
}

# print seqs stored in hash
foreach my $id (keys %hash) {
    $out->write_seq($hash{$id}{'seq'});
}


