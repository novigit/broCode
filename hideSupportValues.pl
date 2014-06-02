#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::TreeIO;

my ($tree_in,$format, $cutoff);
GetOptions('t|tree=s'   => \$tree_in,
	   'f|format=s' => \$format,
	   'c|cutoff=s' => \$cutoff,
);

# objects
my $in = Bio::TreeIO->new(-format => $format, 
			  -file   => $tree_in);
my $out = Bio::TreeIO->new(-format => $format, 
			   -fh     => \*STDOUT);

# loop through trees
while( my $tree = $in->next_tree ) {
    
    # loop through nodes
    for my $node ($tree->get_nodes ) {
	
	# skip terminal nodes (leafs)
	next if ($node->is_Leaf);

	# delete node id (bootstrap value) if lower then 70
	$node->id('') if ($node->id < $cutoff);  
#	print $node->bootstrap, "\n" if ($node->bootstrap < 70);
    }
    
    # write new tree to out object (newick file, STDOUT)
    $out->write_tree($tree);
}

=head1 USAGE

hideSupportValues.pl -t|--tree [input tree] -f|--format [tree format] -c|--cutoff [support value]

=head1 EXAMPLE

hideSupportValues.pl -t|--tree 16S.newick -f|--format newick -c|--cutoff 70

=head1 SYNOPSIS

Hides / Omits support values (bootstraps, posterior probabilities) below the specified value in the provided tree

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)
