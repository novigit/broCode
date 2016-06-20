#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Pod::Usage;

=head1 USAGE

    getGCcontent.pl [fasta]

=head1 SYNOPSIS

    Reports for each sequence in a fasta file, its GC-content
    Also reports the total overall GC-content

=head1 AUTHOR

    Joran Martijn (joran.martijn@icm.uu.se)

=cut

pod2usage (-message => "not enough arguments") unless @ARGV == 1;

my ($fasta) = @ARGV;
my $in  = Bio::SeqIO->new(-format => 'fasta',
			  -file   => $fasta);

my $total_length = 0;
my $total_g_count = 0;
my $total_c_count = 0;

while (my $seqobj = $in->next_seq) {

    # calculate per sequence
    my $id  = $seqobj->id;
    my $seq = $seqobj->seq;
    my $len = $seqobj->length;
    my $g_count = (uc($seq) =~ tr/G/G/);
    my $c_count = (uc($seq) =~ tr/C/C/);
    my $gc = sprintf("%.3f",  ($g_count + $c_count) / $len * 100);

    # report per sequence
    print $id, "\t", $gc, "\n";

    # add to total
    $total_length  += $len;
    $total_g_count += $g_count;
    $total_c_count += $c_count;

}

# calculate overall gc content
my $overallgc = sprintf("%.3f",  ($total_g_count + $total_c_count) / $total_length * 100);
print $fasta, "\t", "Overall:", "\t", $overallgc, "\n";
