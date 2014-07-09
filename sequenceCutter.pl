#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;

=head1 USAGE

sequenceCutter.pl -f|--fasta [.fasta] -i|--id [query] -s|--start [integer] -e|--end [integer] > [out.fasta] 

=head1 SYNOPSIS

Select a contig or sequence in your fasta file, and extract a subsequence, discarding the rest.

=head1 TODO

Now you have to select the subsequence you want to KEEP.
Preferentially want to select the subsequence you want to DISCARD.

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)

=cut


my ($fasta, $id, $start, $end);
GetOptions('f|fasta=s' => \$fasta,
	   'i|id=s'    => \$id,
	   's|start=s' => \$start,
	   'e|end=s'   => \$end);

# in and out objects
my $in  = Bio::SeqIO->new(-file   => $fasta,
			  -format => 'fasta');
my $out = Bio::SeqIO->new(-format => 'fasta',
			  -fh     => \*STDOUT);

# loop over sequences
while ( my $seq = $in->next_seq ) {

    # edit sequence of interest
    if ($seq->id =~ m/$id/) {
	# cut sequence and create object that you can print
	my $subseq = $seq->subseq($start, $end);
	my $subseq_obj = Bio::Seq->new(-seq => $subseq,
				       -id  => $seq->id);
	$out->write_seq($subseq_obj);
    }
    # print all other sequences without edit
    else {
    	$out->write_seq($seq);
    }
}
