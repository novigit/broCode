#!/usr/bin/perl -w
use strict;
use Getopt::Long;

=head1 NAME

countTaxaInCluster.pl

=head1 USAGE

countTaxaInCluster.pl
 
 --amount -n [desired no. of sequences in the COGs] 
 --dir    -d [directory with proteomes, used for the generation of the COGs] 
 --clus   -c [.clus file]

=head1 SYNOPSIS

The aim of this script is to count how many times each taxon is represented in each cluster or COG. It will only scan those COGs which have the desired number of sequences in them. For example, if you specify '-n 24', it will count how many times each taxa is represented in each COG, for all COGs that have a total of 24 sequences. 

First, it creates a hash linking gene name with the genome (taxon) name, using the proteome directory (--dir or -d)

Then, it will parse the .clus file, selecting only those clusters with the desired amount of sequences (--amount or -n). It will then count in the selected clusters, how much a given taxon is present in a hash of hashes:

 Key        Key       Value
 cluster1 - taxon 1 - 2 times
 cluster1 - taxon 2 - 1 time
 cluster1 - taxon 3 - 3 times
 cluster2 - taxon 1 - 1 time
 etc

Lastly it prints the hash of hashes in a table

=head1 COMMENTS

Would like to update to also annotate which clusters are counted in the out file

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)

=cut

die "usage: countTaxaInCluster.pl", "\n", 
"--amount -n\t[desired no. of taxa in COG]", "\n", 
"--dir\t -d\t[directory with cluster]", "\n", 
"--clus\t -c\t[.clus file]", "\n" unless @ARGV == 6;

my @genomes;
my %hash;

my ($amount, $dir, $clusters);
GetOptions('n|amount=s' => \$amount,
	   'd|dir=s'    => \$dir,
	   'c|clus=s'   => \$clusters);

opendir DIR, "$dir";
foreach my $file (readdir DIR) {
    
    unless ($file =~ m/^\./) {

	my ($genome) = $file =~ m/(\S+).faa.fix/;
	push @genomes, $genome;
	# print $genome, "\n";
	# print $file, "\n";
	open FASTA,'<', $dir.'/'.$file;
    	while (<FASTA>) {
    	    if (/^\>/) {
    		my ($gene) = /\>(\S+)/;
    		$hash{$gene} = $genome;
    	    }
    	}
    	close FASTA;
    }
}
closedir DIR;

# foreach my $key (keys %hash) {
#     print $hash{$key}, "\n";
# }

my $clustc = 0;
my %count;

open CLUSTERS, "<$clusters";
while (<CLUSTERS>) {
    my @cluster = split;
    my $size = @cluster;
    if ($size == $amount) {
	$clustc++;
	# print $_;
	foreach my $gene (@cluster) {
            my $genome = $hash{$gene};
	    # print $clustc, "\t", $genome,"\n";
	    $count{$clustc}{$genome}++;
	}
    }
}
close CLUSTERS;

foreach my $genome (@genomes) {
    print $genome, "\t";
}
print "\n";

foreach my $cluster (keys %count) {
    foreach my $genome (@genomes) {
	if ($count{$cluster}{$genome}) {print $count{$cluster}{$genome}}
	else                           {print "0"}
	print "\t";
    }
    print "\n";
}


