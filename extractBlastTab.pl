#!/usr/bin/perl -w

=head1 NAME

    extractBlastTab.pl - extract a subset of a tabular BLAST file based on query names

=head1 USAGE

    extractBlastTab.pl -b|--blasttab <tabular_blast_file> -l|--list <queries.list> > <out.blast>

=head1 DESCRIPTION

    Takes a list with query names as input and searches a tabular blastfile for lines with those queries.
    The script also reports to STDERR the progress of how much queries have been parsed so far

=head1 AUTHOR

    Joran Martijn (joran.martijn@icm.uu.se)

=cut

use strict;
use Getopt::Long;

# options
my ($blasttab, $list);
GetOptions(
    "b|blasttab=s" => \$blasttab,
    "l|list=s"     => \$list,
);

# index the queries in input list
my %hash;
my $qcount = 0;
open LIST, "<", $list;
while (<LIST>) {
    chomp;
    $hash{$_}++;
    $qcount++;
}
close LIST;

# parse tabular blastfile
my $hcount = 0;
my $prev_query = "";
open BLASTTAB, "<", $blasttab;
while (<BLASTTAB>) {

    # capture query name in tabular blast file
    my ($query) = split;
    chomp $query;

    # if query name in blast file matches any in the list
    if (exists $hash{$query}) {

	# report progress to STDERR
	$hcount++ unless ($query eq $prev_query);
	# '%s' to print strings,
	# '%d' to print integer, 
	# '\r' to overwrite instead of newline
	printf STDERR "%s%d%s\r", "Extracted ", $hcount, " out of $qcount queries in $list";

	# print line to out.blast
	print STDOUT $_;
    }

    # store query name to check if next query name has been searched yet
    # this is useful for the counter
    $prev_query = $query;
}
close BLASTTAB;

# Announce finish
print STDERR "Done!"
