#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my ($gbk, $rna) = @ARGV;

my $in  = Bio::SeqIO->new(-format => 'genbank',
			  -file   => $gbk);
my $out = Bio::SeqIO->new(-format => 'fasta',
			  -fh     => \*STDOUT);


# type of rRNA to search
my $query;
if    ($rna eq '16S') { $query = '16S|ssu|small\ subunit' }
elsif ($rna eq '23S') { $query = '23S|lsu|large\ subunit' } 
elsif ($rna eq '5S' ) { $query = '5S'                     }

# loop over records
while (my $seq = $in->next_seq) {

    # get definition or description of sequence
    my $def = $seq->desc;
    next if ($def =~ m/plasmid/i);

    # get species name
    my $strain = $seq->species->node_name;

    # loop over features
    for my $ft ($seq->get_SeqFeatures) {

	# skip all non-rRNA genes
	next unless ($ft->primary_tag eq 'rRNA');

	# get product names
	my ($prod) = $ft->get_tag_values('product') if ($ft->has_tag('product'));
	next unless ($prod =~ m/$query/i);
	# print $prod, "\n";

	# make fasta header
	my $hdr;
	$hdr = $strain . '_notdefined_' . $rna;
	$hdr = $strain . '_complete_'   . $rna if     ($def =~ m/complete/i);
	$hdr = $strain . '_wgs_'        . $rna if     ($def =~ m/whole/i   );
	$hdr =~ s/\s/_/g;

	# get sequence of rRNA genes
	my $ft_seq = $ft->seq;

	# provide $ft_seq with new id and description and print out
    	$ft_seq->id($hdr);
    	$ft_seq->desc('');
    	$out->write_seq($ft_seq);
    }
}
