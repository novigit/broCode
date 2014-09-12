#!/usr/bin/perl

=head1 NAME

    selectSequences.pl - extract or remove sequences from a FASTA file

=head1 USAGE

    selectSequences.pl -i|--fasta [in_fasta] -l|--list [query_list] -q|--queries [query1,query2,query3,etc] -m|--mode [extract|remove] > [out_fasta]

    Use -l or -q
    -q: Query can be a partial match of sequence header
    -l: Query has to be exactly the same as sequence header

=head1 DESCRIPTION

    Extracts a sequence or sequences of interest from a FASTA file.
    The user provides a query or list of queries.
    A query does not have to be an exact match to fasta header, partial works.

=head1 AUTHOR

    Joran Martijn (joran.martijn@icm.uu.se)

=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

pod2usage (-msg => "Not enough arguments") unless @ARGV >= 6;

my ($fasta, $list, $queries, $mode);
GetOptions("i|fasta=s"   => \$fasta,
	   "l|list=s"    => \$list,
	   "q|queries=s" => \$queries,
	   "m|mode=s"    => \$mode);

my $in  = Bio::SeqIO->new(-format => 'fasta',
			  -file   => $fasta);
my $out = Bio::SeqIO->new(-format => 'fasta',
			  -fh     => \*STDOUT);

# If '-q|--queries' was used, create @queries
my @queries;
if ($queries) {
    @queries = split(/,/, $queries);
}
# If '-l|--list'    was used, create %queries
my %queries;
if ($list) {
    open LIST, "<$list";
    while (<LIST>) {
	chomp;
	my $query = $_;
	$queries{$query}++;
    }
    close LIST;
}

# State global $print
my $print;

while (my $seq = $in->next_seq) {
    my $id = $seq->id;
    $id .= $seq->desc;

    # if 'extract' mode, no  sequence  is  printed by default
    $print = 0 if ($mode eq 'extract');
    # if 'remove'  mode, all sequences are printed by default
    $print = 1 if ($mode eq 'remove' );

    # switch print 'on' or 'off' if query matches a sequence header
    print_switch(      \@queries  , $id, $mode) if ($queries);
    print_switch([ keys %queries ], $id, $mode) if ($list)   ;

    # write sequence if print is 'on' 
    $out->write_seq($seq) if $print;
}

sub print_switch {
    my ($var, $id, $mode) = @_;

    # '@$var' either returns '@queries' or 'keys %queries'
    foreach my $query (@$var){

	# turn print 'on'  if query matches sequence header if 'extract' mode
	$print = 1 if ($id =~ m/$query/ && $mode eq 'extract');
	# turn print 'off' if query matches sequence header if 'remove'  mode
	$print = 0 if ($id =~ m/$query/ && $mode eq 'remove' );
    }
    return $print;
}
