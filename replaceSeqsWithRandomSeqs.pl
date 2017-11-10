#!/usr/bin/perl

# some quick documentation:
# this script uses Bio::Factory::EMBOSS module to load and use the 'makeprotseq' program to create random sequence
# the properly load it, make sure that both makeprotseq and another emboss tool 'wossname' are installed in your PATH

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use Bio::Factory::EMBOSS;

my ($fasta, $names);
GetOptions('f|fasta=s' => \$fasta,
	   'n|names=s' => \$names);

# load input fasta
my $in  = Bio::SeqIO->new(-file   => $fasta,
			 -format => 'fasta');
# load output fasta
my $out = Bio::SeqIO->new(-format => 'fasta',
			  -fh     => \*STDOUT);
# load list of sequences you wish to replace
my %seqs_to_replace;
open LIST, "<", $names;
while (<LIST>) {
    	chomp;
	$_ =~ s/^>// if (/^>/);
	$seqs_to_replace{$_}++;
}

# load 'makeprotseq' program
# or, get an EMBOSS application object from the factory
my $f = new Bio::Factory::EMBOSS;
my $makeprotseq = $f->program('makeprotseq');

# loop over sequences
while (my $seq_obj = $in->next_seq) {

    # get sequence id 
    my $id = $seq_obj->id;

    # if sequence is one that you wish to replace
    if ( $seqs_to_replace{$id} ) {

	# get sequence length
    	my $seq_length = $seq_obj->length;

    	# makeprotseq arguments
    	my %flags = ( 
    	    -amount => 1,
    	    -length => $seq_length,
    	    -outseq => 'raw::stdout',
    	    );

    	# run makeprotseq
    	my $random_seq = $makeprotseq->run({%flags});
	$random_seq =~ s{\n}{}g;

    	# catch random sequence in new sequence object
    	my $random_seq_obj = Bio::PrimarySeq->new(-seq        => $random_seq,
						  -display_id => 'random_seq_' . $id,
						  -alphabet   => 'protein');
	
    	# print random sequence
    	$out->write_seq($random_seq_obj);
    }

    # if it isn't, write to outfile unchanged
    else {
    	$out->write_seq($seq_obj)
    }
}
