#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

=head1 USAGE

clusters2matrix.pl [cluster directory] > [clusters.matrix]

=head1 SYNOPSIS

Simply counts for all COGs in a given directory the amount of proteins of a taxon are present. 
Final output:

       taxon1   taxon2    taxon3   ..
COG1     0        0         1      ..
COG2     1        1         0      ..
COG3     1        2         1      ..
COG4     1        1         1      .. 
..      ..       ..        ..      .. 

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se) 

=cut

my ($dir) = @ARGV;

die "usage: clusters2matrix.pl <cluster directory> > <clusters.matrix>\n" unless @ARGV == 1;

my %count;
my %taxa;

opendir DIR, "$dir";
foreach my $faa (readdir DIR) {
    next if ($faa =~ m/^\./);
    print STDERR "Reading $faa ... \n";

    my ($cog) = $faa =~ m/(\S+)\./;

    my $in = Bio::SeqIO->new(-file   => "<$dir$faa",
			     -format => 'fasta');
    while (my $seq = $in->next_seq){
	my $taxon = $seq->id;
	$count{$cog}{$taxon}++;
	$taxa{$taxon}++;
    }
}
closedir DIR;

my @taxa = (keys %taxa);

#print "\t";

print "\t", join "\t", @taxa, "\n";
# foreach my $taxon (@taxa) {
#     print $taxon, "\t"; 
# }
# print "\n";

foreach my $cog (keys %count) {
    print $cog;
    foreach my $taxon (@taxa) {
	if (exists $count{$cog}{$taxon}) {print "\t", $count{$cog}{$taxon}}
	else                             {print "\t", "0"}
    }
    print "\n";
}

