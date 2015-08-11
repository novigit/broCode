#!/usr/bin/perl

=head1 NAME

    selectSequences.pl - extract or remove sequences from a FASTA file

=head1 USAGE

    selectSequences.pl -i|--fasta [in_fasta] -f|--format [fasta/fastq] (default: fasta) -l|--list [query_list] -q|--queries [query1,query2,query3,etc] -m|--mode [extract|remove] -e|--exact -n|--no-description > [out_fasta]

    Use -l or -q

=head1 DESCRIPTION

    Extracts or remove a sequence or sequences of interest from a FASTA file.
    The user provides a query or list of queries.
    A query does not have to be an exact match to fasta header, partial works.
    If the -e mode is used, and exact match is required. This is much faster.

=head1 AUTHOR

    Joran Martijn (joran.martijn@icm.uu.se)
    Anders Lind (anders.lind@icm.uu.se)

=head1 UPDATE

    Anders Lind - Added the possibility to run it on fastq files, 2014-10-02
    Anders Lind - Added the possibility to run exact matches, and to skip looking in the description. 2015-02-06
    Joran Martijn - Added a small fix for when your query contains the ">" sign. The $id does not contain the ">", and thus will not find the query if it contains a ">". 2015-02-20
    Joran Martijn - Fixed issue with concatenating the 'id' and the 'desc', which caused the loss of a space. Added ' ' in concatenation to compensate. 2015-02-23
    Anders lind - Made 'extract' the default mode, and fixed some of the error messages for input flasgs. 2015-06-12
    Anders Lind - Rewrote the matching system and added first hit flag for speed. 2015-06-12

=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

my ($fasta, $list, $queries, $exact,$noDesc,$tophits);
my $format = 'fasta';
my $mode = 'extract'; #Default mode is 'extract'

GetOptions("i|fasta=s"         => \$fasta,
	   "l|list=s"          => \$list,
	   "q|queries=s"       => \$queries,
	   "f|format=s"        => \$format, 
	   "m|mode=s"          => \$mode,
	   "e|exact"           => \$exact,
	   "n|no-description"  => \$noDesc,
	   't|tophits'         => \$tophits,
);

$format = lc($format); #change to lower-case if someone specified with uppercase
pod2usage (-msg => "No input fasta specified") unless $fasta;
pod2usage (-msg => "No list OR query specified") unless ($list || $queries);
pod2usage (-msg => "Format can only be 'fasta' or 'fastq'") unless ($format eq "fastq" || $format eq "fasta");

my $in  = Bio::SeqIO->new(-format => $format,
			  -file   => $fasta);
my $out = Bio::SeqIO->new(-format => $format,
			  -fh     => \*STDOUT);

my %queries;
# If '-q|--queries' was used, create @queries
my @queries;
if ($queries) {
    @queries = split(/,/, $queries);
    foreach (@queries) {
	$_ =~ s/^>// if (/^>/);
	$queries{$_}++;
    }
    
}
# If '-l|--list'    was used, create %queries
if ($list) {
    open LIST, "<$list";
    while (<LIST>) {
	chomp;
	$_ =~ s/^>// if (/^>/);
	$queries{$_}++;
    }
    close LIST;
}

# State global $print
my $print;
my $counter = 0;
while (my $seq = $in->next_seq) {
    my $id = $seq->id;
    if ($seq->desc) {
	$id .= ' ' . $seq->desc unless ($noDesc); #Skip looking in description
    }
    # if 'extract' mode, no  sequence  is  printed by default
    $print = 0 if ($mode eq 'extract');
    # if 'remove'  mode, all sequences are printed by default
    $print = 1 if ($mode eq 'remove' );

    if ($exact) {
	#print STDERR "Exact in use\n";
	$print = 1 if ($queries{$id} && $mode eq 'extract');
	$print = 0 if ($queries{$id} && $mode eq 'remove');
    }
    else {
    # switch print 'on' or 'off' if query matches a sequence header
	print_switch([ keys %queries ], $id, $mode);
    }
    # write sequence if print is 'on' 
    $out->write_seq($seq) && $counter++ if $print;
    last if ( ($tophits) && ($counter >= scalar(keys %queries)) ); #Stop looking when all found
    
}

sub print_switch {
    my ($var, $id, $mode) = @_;
    foreach my $query (@$var){
	# turn print 'on'  if query matches sequence header if 'extract' mode
	$print = 1 if ($id =~ m/\Q$query\E/ && $mode eq 'extract');
	# turn print 'off' if query matches sequence header if 'remove'  mode
	$print = 0 if ($id =~ m/\Q$query\E/ && $mode eq 'remove' );
    }
    return $print;
}
