package Phylogeny::Tree::Node;
use strict;
use warnings;
use feature ':5.10';

use Moose;

has _children => (
    is => 'ro',
    #isa => 'ArrayRef',
    traits => ['Array'],
    handles => {
        children => 'elements',
        add_child => 'push',
    },
    default => sub { [] },
);

has ancestor => (
    is => 'ro',
    required => '',
);

has name => (
    is => 'ro',
    #isa => 'Str',
    required => '',
);

sub is_leaf {
    my $self = shift;
    return $self->children == 0;
}

sub get_all_children {
    my $self = shift;
    return $self if $self->is_leaf;
    return map { $_->get_all_children } $self->children;
}

sub get_all_nodes {
    my $self = shift;
    return $self if $self->is_leaf;
    return ($self, map { $_->get_all_nodes } $self->children);
}

sub get_all_internal_nodes {
    my $self = shift;
    return if $self->is_leaf;
    return ($self, map { $_->get_all_internal_nodes } $self->children );
}

sub get_splits {
    my $self = shift;
    my @all_nodes = @{ $self->{cnames} };

    my @splits;
    for my $node ( $self->get_all_internal_nodes ) {
        my @below = @{ $node->{cnames} };
        my %below = map { $_,1 } @below;
        my @above = grep { ! exists $below{$_} } @all_nodes;

        next unless @below > 1 && @above > 1;

        push @splits, Phylogeny::Tree::Split->new( \@above, \@below );
    }
    return @splits;
}

__PACKAGE__->meta->make_immutable;

1;
