#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my ($gbk, $cds) = @ARGV;

my $in  = Bio::SeqIO->new(-format => 'genbank',
			  -file   => $gbk);
my $out = Bio::SeqIO->new(-format => 'fasta',
			  -fh     => \*STDOUT);


# loop over records
while (my $seq = $in->next_seq) {

    # get definition or description of sequence
    my $def = $seq->desc;

    # get species name
    my $strain = $seq->species->node_name;

    # loop over features
    for my $ft ($seq->get_SeqFeatures) {

	# skip all non-CDS features
	next unless ($ft->primary_tag eq 'CDS');

	# skip CDS features that have pseudo tags
	next if ($ft->has_tag('pseudo'));

	# get CDS accession
	my ($acc) = $ft->get_tag_values('protein_id') if ($ft->has_tag('protein_id'));
	# print $acc, "\n";

	# get CDS product names
	my ($prod) = $ft->get_tag_values('product') if ($ft->has_tag('product'));
	# print $prod, "\n";

	# make fasta header
	my $hdr;
	$hdr = $strain . '_i_' . $acc . '_l_' . $prod;
	$hdr =~ s/\s/_/g;
	# print $hdr, "\n";

	# get sequence of CDS
	my ($aa_seq) = $ft->get_tag_values('translation') if ($ft->has_tag('translation'));
	# print $aa_seq, "\n";

	# provide $ft_seq with new id and description and print out
	my $outseq_obj = Bio::Seq->new(-seq  => $aa_seq,
				       -id   => $hdr,
				       -desc => '',
	);

	# my $ft_seq;
	# $ft_seq->seq($aa_seq);
    	# # $ft_seq->id($hdr);
    	# # $ft_seq->desc('');
    	$out->write_seq($outseq_obj);
    }
}
