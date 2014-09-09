#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Pod::Usage;

=head1 USAGE

    getGCcontent.pl [fasta]

=head1 SYNOPSIS

    Reports for each sequence in a fasta file, its GC-content

=head1 AUTHOR

    Joran Martijn (joran.martijn@icm.uu.se)

=cut

pod2usage (-message => "not enough arguments") unless @ARGV == 1;

my ($fasta) = @ARGV;
my $in  = Bio::SeqIO->new(-format => 'fasta',
			  -file   => $fasta);

while (my $seqobj = $in->next_seq) {
    my $id  = $seqobj->id;
    my $seq = $seqobj->seq;
    my $len = $seqobj->length;
    my $g_count = (uc($seq) =~ tr/G/G/);
    my $c_count = (uc($seq) =~ tr/C/C/);
    my $gc = sprintf("%.3f",  ($g_count + $c_count) / $len * 100);

    print $id, "\t", $gc, "\n";

}
