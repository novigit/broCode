#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my ($gbk) = @ARGV;

my $in = Bio::SeqIO->new(-file   => $gbk,
			 -format => 'genbank');

# genbank record categories
my $wgs='whole genome shotgun';
my $complete='complete genome';
my $plasmid='plasmid';
my $chrmsm='complete sequence';

# array that will store all accessions
my @accessions;

while (my $seq = $in->next_seq) {
    my $def = $seq->desc;

    # skip plasmids
    next if ($def =~ m/$plasmid/i);

    # get accession from complete genomes
    if ($def =~ m/$complete|$chrmsm/i) {
    	my $acc = $seq->accession_number;
	push @accessions, $acc;
    }
    
    # get accessions from whole genome shotgun projects
    elsif ($def =~ m/$wgs/i) {

	# get Genbank fields (like COMMENT, REFERENCE, DBLINK, WGS, WGS_SCAFLD etc)
	# WGS_SCAFLD is for RefSeq, WGS is for WGS
	my $anns = $seq->annotation;

	# get annotation of WGS_SCAFLD field(s) and WGS field(s)
	my @wgs_anns = $anns->get_Annotations('wgs'); 
	my @wgs_scafld_anns = $anns->get_Annotations('wgs_scafld');
	my @anns = (@wgs_anns, @wgs_scafld_anns);
	
	# retrieve accessions
	for my $val (@anns) {
	    my $wgs_scafld = $val->display_text;
	    
	    # if multiple scaffolds
	    if ($wgs_scafld =~ m/-/) {
		my ($acc_start, $acc_end) = split(/-/, $wgs_scafld);
		
		# get the full range of accession numbers of scaffolds
		my ($id)    = $acc_start =~ m/(^\D+)/;
		my ($start) = $acc_start =~ m/(\d+)$/;
		my ($end)   = $acc_end   =~ m/(\d+)$/;
		push @accessions, map{ $id.$_ }($start..$end);
	    }

	    # if single scaffold
	    else { 
		push @accessions, $wgs_scafld 
	    }
	}
    }
}

#print all stored accessions
print join("\n", @accessions), "\n";
