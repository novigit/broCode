#!/usr/bin/perl -w
use strict;
use Bio::NEXUS;
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;
=head1 NAME

    newick2nexus.pl

=head1 USAGE

    newick2nexus.pl -i|--infile [.newick] > [.nexus]
    
    -i|--infile  COMPULSORY   The newick-file you want to convert
    -r           OPTIONAL     Convert '-' into '_' in your newick-file

    Assumes taxon names do not include spaces

=head1 AUTHOR

    Lionel Guy, lionel.guy@icm.uu.se (original)
    Joran Martijn, joran.martijn@icm.uu.se (edit)

=cut

## Read from file
my $infile;
my $remove;
my ($fh, $filename) = tempfile();

GetOptions(
    'i|infile=s' => \$infile,
    'r|remove-non-alphanumeric-characters' => \$remove,
);
my $INCON = \*STDIN;
if ($infile){
    open $INCON, '<', $infile or die;
}

my ($tree_str) = <$INCON>;
chomp $tree_str;

## Remove dashes. Hopefully prevents the nexus module to screw up
$tree_str =~ s/-/_/g if $remove;

## Create an empty Trees Block, and then add a tree to it
my $trees_block   = new Bio::NEXUS::TreesBlock('trees');

$trees_block->add_tree_from_newick($tree_str, "my_tree");

## Create new Bio::NEXUS object
my $nexus_obj = new Bio::NEXUS;
$nexus_obj->add_block($trees_block);
$nexus_obj->write($filename);

## Fine tune and print
open NEX, '<', $filename;
while (<NEX>){
    chomp;
    if (/(\s+)(TAXLABELS\s+.*)$/){
    	my $spaces = $1;
    	my @taxlabels = split(/\s+/, $2);
#	if ($taxlabels[-1] =~ /;$/){
	    $taxlabels[-1] =~ s/\;$/\n;/;
#	}
    	print $spaces . join("\n$spaces", @taxlabels) . "\n";
    }
    elsif (/\s+TREE my_tree =/){
    	s/inode//g;
    	s/\';\'//g;
    	print "$_\n";
    }
    else {
    	print "$_\n";
    }
}
