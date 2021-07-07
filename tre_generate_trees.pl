#!/usr/bin/perl
package tre_generate_trees;
use strict;
use warnings;
use feature ':5.10';
use Moose;
use MooseX::Types::Path::Class;
with 'MooseX::Getopt';

use Data::Printer;

sub run {
    my $self = shift;

    my @taxa = qw/ A B C D /;

    my @trees = $self->make_trees( @taxa );
    for my $tree ( @trees ) {
        say $self->tree_string( $tree );
    }
}

sub tree_string {
    my $self = shift;
    my $tree = shift;
    if ( ! ref $tree ) {
        return $tree;
    }
    my $str = join ',', map { $self->tree_string($_) } @$tree;
    return "($str)"
}

sub make_trees {
    my $self = shift;
    my @elems = @_;

    if ( @elems == 2 ) {
        return ([@elems]);
    }

    my @all = @elems;

    my @trees;
    for my $idx (0..$#elems) {
        @elems = @all;
        my $e = splice @elems, $idx, 1;
        my @subtrees = $self->make_trees(@elems);
        push @trees, map { [$e, $_] } @subtrees;
    }
    return @trees;
}

tre_generate_trees->new_with_options->run unless caller;

1;

