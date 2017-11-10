#!/usr/bin/perl -w
use strict;

=head1 NAME

theSubstitutor.pl - Searches and replaces a text file according to a mapping file

=head1 USAGE

theSubstitutor.pl -i <input> -m <mapping_file> > <output>

Make sure the mapping file is like this
<search_pattern1>\t<replace_pattern1>
<search_pattern2>\t<replace_pattern2>

=head1 DESCRIPTION

Stores the mapping file in a hash and searches per line for your search patterns

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)

=cut

# libraries
use Getopt::Long;
use Pod::Usage;

# state options
my ($map, $input);
GetOptions(
    'i|input=s' => \$input,
    'm|map=s'   => \$map,
    );
pod2usage (-msg => "Not enough arguments") unless $input && $map;

# datastructures
my %mapping;

# load mapping into memory
open MAP, '<', $map;
while (my $line = <MAP>){
    chomp $line;
    my ($search, $replace) = split '\t', $line;
    # print $search, "\t", $replace, "\n";
    $mapping{$search}=$replace;
}
close MAP;

my $regex = join "|", keys %mapping; # '|' for the OR in REGEX
#print $regex, "\n";

$regex = qr/$regex/; # qr/STRING/ returns STRING interpreted as REGEX pattern
#print $regex, "\n";

# replace text in input
open INPUT, '<', $input;
while (my $line = <INPUT>){
    if (my ($hit) = $line =~ m/($regex)/) {
	#print $hit, "\n"; 
	#print $hit, "\t", $mapping{$hit}, "\n";
	$line =~ s/$hit/$mapping{$hit}/g;
	print $line;
    }
    else {
    	print $line;
    }
}
close INPUT;

