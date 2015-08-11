#!/usr/bin/perl -w
use strict;
use List::MoreUtils qw(any);

my ($rnt) = @ARGV;

my $ssu_start;
my $lsu_end;
my $ssu_end;
my $lsu_start;
my $distance;
my $lsu_ssu_set_count = 0;

my @distances;

open RNT, "<$rnt";
while (<RNT>) {
    # skip first 3 lines
    next unless ($. > 3);
    chomp;

    # get coordinates, strand info and product names
    my @line = split "\t", $_;
    my ($coord, $strand, $prod) = ($line[0], $line[1], $line[8]);
    # print $coord, $strand, $prod, "\n";

    my ($start,$end) = split /\.\./, $coord;
    
    if ($strand eq '-' && $prod =~ m/23S/) {
	$lsu_end = $end;
	# print $lsu_end, "\n";
    }

    if ($strand eq '-' && $prod =~ m/16S/) {
	$ssu_start = $start;
	# print $ssu_start, "\n";
	if ($lsu_end) {
	    $distance = $ssu_start - $lsu_end;
	    $lsu_ssu_set_count++ if ($distance < 1000);
	    # print $distance, "\n";
	    push @distances, $distance;
	}
	else { print "LSU before SSU on negative strand\n" }
    }

    if ($strand eq '+' && $prod =~ m/16S/) {
	$ssu_end = $end;
	# print $ssu_end, "\n";
    }

    if ($strand eq '+' && $prod =~ m/23S/) {
	$lsu_start = $start;
	# print $lsu_start, "\n";
	if ($ssu_end) {
	    $distance = $lsu_start - $ssu_end;
	    $lsu_ssu_set_count++ if ($distance < 1000);
	    # print $distance, "\n";
	    push @distances, $distance;
	}
	else { print "SSU before LSU on negative strand\n" }

   }
}
close RNT;

# print join("\n",@distances), "\n";
# print $rnt, "\n";
print "$rnt has neighboring SSU and LSU genes\n" if (any { $_ < 1000 } @distances);
# print "$lsu_ssu_set_count\n";
# print "#####################\n";
