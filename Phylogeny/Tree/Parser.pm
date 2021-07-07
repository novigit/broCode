package Phylogeny::Tree::Parser;
use strict;
use warnings;
use feature ':5.10';
no if $] >= 5.018, warnings => "experimental::smartmatch";

use Carp;

use Phylogeny::Tree::Node;

sub parse {
    my $self = shift;
    my $str = shift;

    # Removing branchlenghts and bootstraps
    $str =~ s/:\d+\.?\d*//g;
    $str =~ s/\)\K\d+//g;

    # Removing genenames
    #$str =~ s/\|[0-9A-F]+//ig;
    #$str =~ s/_\w+//g;

    my $tokens = $self->tokenize( $str );
    my $root   = $self->get_tree( $tokens );

    return $self->preprocess_tree( $root );
}

sub preprocess_tree {
    my $self = shift;
    my $tree = shift;

    for my $child ( $tree->get_all_children ) {
        my $name = $child->name;
        while ( $child ) {
            push @{ $child->{cnames} }, $name;
            $child = $child->ancestor;
        }
    }

    return $tree;
}

sub tokenize {
    my $self = shift;
    my $str = shift;
    # Tokenizer
    my @tokens;
    my $curname = '';
    for ( split //, $str ) {
        when ( '(' ) { push @tokens, { type => 'BRANCH_START' } }
        when ( ')' ) {
            if ( $curname ) {
                push @tokens, { type => 'NODE', name => $curname };
                $curname = '';
            }
            push @tokens, { type => 'BRANCH_END'   } }
        when ( ',' ) {
            if ( $curname ) {
                push @tokens, { type => 'NODE', name => $curname };
                $curname = '';
            }
        }
        default { $curname .= $_ }
    }

    return \@tokens;
}

sub get_tree {
    my $self = shift;
    my $tokens = shift;

    my $root = Phylogeny::Tree::Node->new;
    my $current_node = $root;
    for my $token ( @$tokens ) {
        given ( $token->{type} ) {
            when ('NODE') {
                my $node = new_node( $current_node, $token->{name} );
            }
            when ( 'BRANCH_START' ) {
                $current_node = new_node( $current_node );
            }
            when ( 'BRANCH_END' ) {
                $current_node = $current_node->{ancestor}
            }
        }
    }

    if ( $root->children == 1 ) {
        ($root) = $root->children;
        $root->{ancestor} = undef;
    }
    return $root;
}

sub new_node {
    my $ancestor = shift;
    my $name     = shift;
    my $node = { ancestor => $ancestor };
    if ( $name ) {
        $node->{name} = $name;
    }
    $node = Phylogeny::Tree::Node->new( $node );
    $ancestor->add_child( $node );
    return $node;
}

1;

