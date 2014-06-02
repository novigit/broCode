#!/usr/bin/perl -w
use strict;

# Parses the bitwise flag field of the sam file and counts happy paired reads, i.e. mapped within expected insert size and correct orientation#

die "usage: parseSAMHappyPairs.pl <SAM file>\n" unless @ARGV == 1;

my $total_count;
my $happy_count;
my $unhappy_count;

open IN, "< $ARGV[0]";

while (<IN>) {
    unless ($_ =~ /^@/) {
	$total_count++;                
	my @line = split '\t', $_;     
	my $flag = $line[1];           
#	print "$flag\n";
	if ($flag == 99 || $flag == 147 || $flag == 83 || $flag == 163) {  
	    # print "Happy Pair\t";
	    # print "$_";
	    $happy_count++;
	}
	else {
	    # print "Unhappy Pair\t";
	    # print "$_";
	    $unhappy_count++;
	}
    }
}
my $happy_pair_count = $happy_count / 2;
my $unhappy_pair_count = $unhappy_count / 2;
my $total_pair_count = $total_count / 2;
my $happy_percentage = ($happy_pair_count / $total_pair_count) * 100;

close IN;

print "Total amount of Pairs: $total_pair_count\n";
print "Amount of Happy Pairs: $happy_pair_count\n";
print "Amount of Unhappy Pairs: $unhappy_pair_count\n";
print "% Happy Pairs: "; printf "%.1f\n", "$happy_percentage";
