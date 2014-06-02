#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;

=head1 NAME

clusters2faa.pl

=head1 SYNOPSIS

Builds clusters (protein fasta files) of orthologous groups from the faMCL normalCluster.clus file. This script will build áll the clusters given in the .clus file, not just the panorthologs.

=head1 COMMENTS

Heavily inspired on Lionel's clustersFaa2core.pl.

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)

=cut

# Usage
die "usage: clusters2faa.pl", "\n",
    "-c|--clus", "\t", "[.clus file]", "\n",
    "-p|--proteomes", "\t", "[directory from which to extract protein sequence]", "\n",
    "-o|--outdir", "\t", "[where to put the cogs]", "\n",
    "-n|--name", "\t", "[the name of your cogs]", "\n" unless @ARGV == 8;

# Options
my ($clus, $prot_dir, $out_dir, $cog_name);
GetOptions('c|clus=s'      => \$clus,
	   'p|proteomes=s' => \$prot_dir,
	   'o|outdir=s'    => \$out_dir,
	   'n|name=s'      => \$cog_name);

# Main data structures
my %genes;

# Parse protein files
opendir DIR, "$prot_dir";
foreach my $fasta (readdir DIR) {
    next if ($fasta =~ m/^\./);
    print STDERR "Reading $fasta ...\n";

    # name species from fasta file
    my ($org) = $fasta =~ m/(\S+).faa.fix/;
    
    # store sequence in two-dimensional hash
    my $fasta_in = Bio::SeqIO->new(-file   => "<$prot_dir$fasta",
				   -format => 'fasta');
    while (my $seq = $fasta_in->next_seq){
	my $id = $seq->id();
	$genes{$id}{'org'} = $org;
	$genes{$id}{'seq'} = $seq;
    }
	
}
closedir DIR;

print "==============LOADED ALL GENES=================\n";

# Parse .clus file
open CLUS, '<', $clus;
print STDERR "Reading $clus ...\n";
my $count = 0;
while (<CLUS>) {  
    chomp;

    # specifying the name of the clusters
    $count++;
    my $name = "$cog_name" . sprintf("%05d", $count);
    
    # specifying outfile
    my $fasta_out = Bio::SeqIO->new(-file   => ">$out_dir$name.faa",
    				    -format => 'fasta');
	
    # print sequences of this cluster out to new fasta file
    my @genes = split "\t", $_;
    foreach my $gene (@genes) {
	
    	# retrieve taxon name and sequence from %genes
    	my $org = $genes{$gene}{'org'};
    	my $seq = $genes{$gene}{'seq'};
    	die "$gene was not found in provided protein files\n" unless $org;
	
    	# print out
    	$seq->id($org);
    	$seq->desc($gene);
    	$fasta_out->write_seq($seq);
	
    }
}
close CLUS;
