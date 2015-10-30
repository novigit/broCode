#!/usr/bin/perl

=head1 NAME

    gff2fasta.pl - converts gff to fasta

=head1 USAGE

    gff2fasta.pl -g|--gff <input.gff> -f|--fasta <input.fasta> > <output.fasta>

=head1 DESCRIPTION

    Takes a gff file (for example barrnap output) and an input fasta file (for example contigs) 
    and returns fasta with sequence of features in input gff
    gff version 3

=head1 AUTHOR

    Joran Martijn (joran.martijn@icm.uu.se)
    
=cut

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::FeatureIO;
use Pod::Usage;

# get options
my ($gff,$fasta);
GetOptions("g|gff=s"   => \$gff,
	   "f|fasta=s" => \$fasta,
);
pod2usage (-msg => "Not enough arguments") unless $gff && $fasta;

# parse gff 
my %gff;
open GFF, '<', $gff;
while(<GFF>) {
    my @ft = split;
    my ($g_seqid,$g_start,$g_end,$g_attributes) = ($ft[0],$ft[3],$ft[4],$ft[8]);

    # make 2dim hash of referenced arrays
    # 2dim hash because several features have the same seqid
    $gff{$g_seqid}{$.} = [$g_start, $g_end, $g_attributes]; 
}
close GFF;

# in_fasta object
my $in_fasta  = Bio::SeqIO->new(-format => 'fasta',
			       -file   => $fasta
    );

# out_fasta object
my $out_fasta = Bio::SeqIO->new(-format => 'fasta',
				-fh     => \*STDOUT
    );

# parse fasta
while (my $seq = $in_fasta->next_seq) {

    my $f_seqid = $seq->id;

    # skip seqs that are not stated in gff
    next unless $gff{$f_seqid};
    # print $f_seqid, "\t";
    
    foreach my $key (keys $gff{$f_seqid}) {
	# print $key, "\t";
	
	## cut in_fasta sequence according to gff coordinates
	# first access hash and dereference referenced array
	my ($f_start,$f_end,$f_attributes) = @{$gff{$f_seqid}{$key}};
	# print $f_seqid, "\t", $f_start, "\t", $f_end, "\n";

	# then extract sequence
	my $subseq = $seq->subseq($f_start, $f_end);

	# and provide attribute as description in new seq object
	my $subseq_obj = Bio::Seq->new(-seq    => $subseq,
				       -id     => $f_seqid,
				       -desc   => $f_attributes,
	);
	$out_fasta->write_seq($subseq_obj);
    }
}
