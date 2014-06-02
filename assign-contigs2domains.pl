#!/usr/bin/perl -w
use strict;
use Getopt::Long;

=head1 NAME

assign-contigs2domains.pl

=head1 USAGE

assign-contigs2domains.pl --domain [ar|bac|euk|vir|nh] --filter [only|majority] --cutoff [0<x<1] [IN.tsv] > [OUT.tsv]

IN:  tab seperated file that lists for each contig how many genes belong to which taxonomic domain and phylym
OUT: tab seperated file with a sub selection of contigs which correspond to a certain domain

=head1 SYNOPSIS

This script assigns contigs to a certain taxonomic domain, based on the fraction of genes which belong to that domain. For example, if the cutoff is set to 0.67, and more than two-thirds of the genes are archaeal, that contig will be assigned as Archaeal. Domain assigned contigs are printed out in the [OUT.tsv] file.

=head1 OPTIONS

--domain => Selects taxonomic domain you want to assign contigs with
--filter => Stringency of assignment. 
    'only'     will assign contigs if they have genes to are assigned to one domain only
    'majority' will assign contigs if the fraction [genes of a selected domain / total no. of genes] exceeds the cutoff
--cutoff => Sets the cutoff the fraction will have to exceed in order to assign a contig to the selected domain

=head1 COMMENTS

Currently this script only works with the output of assign-genes2domains.pl as input

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)

=cut


die "usage: assign-contigs2domains.pl --domain [ar|bac|euk|vir|nh] --filter [only|majority] --cutoff [0<x<1] [domain_count.tsv]\n" unless (@ARGV == 5 || @ARGV == 7);

my ($domain, $filter,$cutoff);
GetOptions("d|domain=s" => \$domain,
           "f|filter=s" => \$filter,
           "c|cutoff=s" => \$cutoff); 

open IN, "<$ARGV[0]";

print scalar <IN>;

while (<IN>) {
    chomp;
    my ($node, @taxa) = split '\t', $_;

    my $sum = 0;
    my $limit = @taxa;
    for (my $i=0; $i<$limit; $i++) {
	# print $taxa[$i], "\n";
	$sum+=$taxa[$i];
    }

    my $nh = pop @taxa;

    if ($filter eq 'only') {
	print $_, "\n" if ($domain eq 'ar'  && $taxa[0]+$taxa[1] == $sum);
	print $_, "\n" if ($domain eq 'bac' && $taxa[2]+$taxa[3] == $sum);
	print $_, "\n" if ($domain eq 'euk' && $taxa[4]          == $sum);
	print $_, "\n" if ($domain eq 'vir' && $taxa[5]          == $sum);
	print $_, "\n" if ($domain eq 'nh'  && $nh               == $sum);
    }

    if ($filter eq 'majority') {
	print $_, "\n" if ($domain eq 'ar'  && $taxa[0]+$taxa[1] > $cutoff * $sum);
	print $_, "\n" if ($domain eq 'bac' && $taxa[2]+$taxa[3] > $cutoff * $sum);
	print $_, "\n" if ($domain eq 'euk' && $taxa[4]          > $cutoff * $sum);
	print $_, "\n" if ($domain eq 'vir' && $taxa[5]          > $cutoff * $sum);
	print $_, "\n" if ($domain eq 'nh'  && $nh               > $cutoff * $sum);
    }
}

close IN;
    

