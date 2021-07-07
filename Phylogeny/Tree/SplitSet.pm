package Phylogeny::Tree::SplitSet;
use strict;
use warnings;
use feature ':5.10';
use Moose;

use Carp;
use Phylogeny::Tree::Split;

# The _splits datastructure is a hashref of arrayrefs where the key is the
# number of left elements in a split and the arrayref is a list of hashrefs
# with a 'split' key for the split and a 'count' key showing the number of
# times for that split.
#
# Example:
#  _splits = {
#      2 => [
#               { split => Phylogeny::Tree::Split(left=>['A,B'],right=>['C,D']), count => 2 },
#               { split => Phylogeny::Tree::Split(left=>['A,C'],right=>['B,D']), count => 10 },
#           ],
#      3 => ....
#  }

has _splits => (
    is       => 'ro',
    #isa      => 'HashRef',
    required => '',
    default => sub { {} },);

has _all_splits => (
    is => 'ro',
    lazy_build => 1,
    clearer => '_invalidate_all_splits',
);

has _species => (
    is => 'ro',
    #isa => 'ArrayRef',
    lazy_build => 1,
    traits => ['Array'],
    handles => {
        species => 'elements',
    }
);

sub _build__species {
    my $self = shift;
    my ($asplit) = $self->splits;
    return [ $asplit->{split}->species ];
}

sub _build__all_splits {
    my $self = shift;
    return [ sort { $a->{split}->order($b->{split}) } map { @$_ } values %{ $self->_splits } ];
}

sub from_string {
    my $self = shift;
    my $str  = shift;

    my $_splits;
    for my $split_str ( split /\n+/, $str ) {
        my ($count, $split) = split /\s/,$split_str,2;
        $split = Phylogeny::Tree::Split->new( $split );
        my $num_left = $split->num_left;
        push @{ $_splits->{$num_left} }, { count => $count, split => $split };
    }

    return Phylogeny::Tree::SplitSet->new( _splits => $_splits );
}

sub stringify {
    my $self = shift;
    my $str = join "\n",
        map { $_->{count} . ' ' . $_->{split}->stringify }
        $self->splits;
    return $str;
}

sub add_split {
    my $self = shift;
    my $new_split = shift;
    my $num       = shift // 1;

    my $num_left = $new_split->num_left;

    my $splits = $self->_splits;
    if ( exists $splits->{$num_left} ) {
        for my $split ( @{ $splits->{$num_left} } ) {
            if ( $new_split->is_identical( $split->{split} ) ) {
                $split->{count} += $num;
                return;
            }
        }
    }

    $self->_invalidate_all_splits;
    push @{ $splits->{$num_left} }, { split => $new_split, count => $num };
}

sub add_tree {
    my $self = shift;
    my $tree = shift;

    my @leaves = map { $_->id } $tree->get_leaf_nodes;
    my @internal_nodes = grep { !$_->is_Leaf } $tree->get_nodes;

    my @splits;
    for my $node ( @internal_nodes ) {
        my @below = map { $_->id } grep { $_->is_Leaf } $node->get_all_Descendents;
        my %below = map  { $_ => 1 } @below;
        my @above = grep { !$below{$_} } @leaves;

        next unless @above && @below;

        my $split = Phylogeny::Tree::Split->new( \@above, \@below );
        $self->add_split( $split );
    }
}

sub prune {
    my $self = shift;
    my $cutoff = shift;

    if ( ! $cutoff || ! defined $cutoff ) {
        croak "Need cutoff value for prune";
    }

    my $_splits = $self->_splits;
    for my $subnum ( keys %$_splits ) {
        my @new_split = grep { $_->{count} >= $cutoff } @{ $_splits->{$subnum} };
        if ( @new_split ) {
            $_splits->{$subnum} = \@new_split;
        }
        else {
            delete $_splits->{$subnum};
        }
    }
    $self->_invalidate_all_splits;
}

sub splits {
    my $self = shift;
    my $cutoff = shift // 0;

    if ( $cutoff ) {
        return grep { $_->{count} >= $cutoff } @{ $self->_all_splits };
    }
    return @{ $self->_all_splits };
}

sub num_splits {
    my $self = shift;
    my $cutoff = shift // 0;

    my $count = grep { $_->{count} >= $cutoff } @{ $self->_all_splits };
    return $count;
}

sub remove_species {
    my $self = shift;
    my @species = @_;

    my $new_split_set = Phylogeny::Tree::SplitSet->new();
    for my $split ( $self->splits ) {
        my $new_split;
        eval {
            $new_split = $split->{split}->remove_species( @species );
        };
        next if $@;
        $new_split_set->add_split( $new_split, $split->{count} );
    }

    return $new_split_set;
}

sub conflict_score {
    my $self   = shift;
    my $other  = shift;
    my $cutoff = shift // 0;
    my $same_species = shift // '';

    if ( ! $same_species ) {
        my %all_species = map { ($_,1) } $self->species;
        $all_species{$_}++ for $other->species;
        my @singles = grep { $all_species{$_} == 1 } keys %all_species;
        if ( @singles ) {
            #printf STDERR "Will remove @singles\n";
            $self  = $self->remove_species( @singles );
            $other = $other->remove_species( @singles );
        }
    }

    my @self_splits  = map { $_->{split} } $self->splits( $cutoff );
    my @other_splits = map { $_->{split} } $other->splits( $cutoff );

    my $num_conflicts = 0;
    for my $self_split ( @self_splits ) {
        for my $other_split ( @other_splits ) {
            if ( ! $self_split->is_compatible( $other_split ) ) {
                $num_conflicts++;
            }
        }
    }
    return $num_conflicts/(@self_splits * @other_splits);
}

__PACKAGE__->meta->make_immutable;

1;
