package Phylogeny::Tree::Split;
no if $] >= 5.018, warnings => "experimental::smartmatch";
use strict;
use warnings;
use feature ':5.10';
use Moose;
use Carp;

has left => (
    is => 'ro',
    #isa => 'ArrayRef[Str]',
    required => 1,
    traits => ['Array'],
    handles => {
        _lefts => 'elements',
        num_left => 'count',
    },
);

has right => (
    is => 'ro',
    #isa => 'ArrayRef[Str]',
    required => 1,
    traits => ['Array'],
    handles => {
        _rights => 'elements',
        num_right => 'count',
    },
);

has stringify => (
    is => 'ro',
    isa => 'Str',
    lazy_build => 1,
);

has _species => (
    is => 'ro',
    isa => 'HashRef',
    lazy_build => 1,
    traits => ['Hash'],
    handles => {
        __species => 'keys',
        has_species => 'exists',
    }
);

# Make sure that the split is in a canonical form. Always the shortest part of
# the split to the left. If they are equally long, the first one in
# lexiographic order is to the left.
around BUILDARGS => sub {
    my $orig = shift;
    my $class = shift;

    my ($left,$right);

    # Extract the left and right components
    if ( @_ == 1 && ref $_[0] eq 'HASH' ) {
        $left = $_[0]{left};
        $right = $_[0]{right};
    }
    elsif ( @_ == 4 ) {
        my %hash = @_;
        ($left, $right) = @hash{qw/ left right /};
    }
    elsif ( @_ == 2 ) {
        ($left, $right) = @_;
    }
    elsif ( @_ == 1 && ! ref $_[0] ) {
        ($left,$right) = map { [split ',',$_] } split ' ', $_[0];
    }
    else {
        return $class->$orig();
    }

    # Sort them
    my @left  = sort @$left;
    my @right = sort @$right;

    if ( @left < 2 || @right < 2 ) {
        croak "Both left and right need to contain at least 2 elements";
    }

    # Build the object
    if ( @left < @right ) {
        return $class->$orig( left => \@left,  right => \@right);
    }
    elsif ( @left > @right ) {
        return $class->$orig( left => \@right, right => \@left );
    }
    elsif ( $left[0] lt $right[0] ) {
        return $class->$orig( left => \@left,  right => \@right);
    }
    elsif ( $left[0] gt $right[0] ) {
        return $class->$orig( left => \@right, right => \@left );
    }

    # Something went wrong
    return $class->$orig();
};

sub _build__species {
    my $self = shift;
    return { map { ($_,1) } ($self->_lefts, $self->_rights) };
}

sub _build_stringify {
    my $self = shift;
    return join ' ', $self->left_string, $self->right_string;
}

sub species {
    my $self = shift;
    return sort $self->__species;
}

sub remove_species {
    no if $] >= 5.018, warnings => "experimental::smartmatch";

    my $self = shift;
    my @species = @_;

    my @left  = grep { ! ($_ ~~ @species) } $self->_lefts;
    my @right = grep { ! ($_ ~~ @species) } $self->_rights;

    return Phylogeny::Tree::Split->new(\@left, \@right);
}

sub left_string {
    my $self = shift;
    return join ',', $self->_lefts;
}

sub right_string {
    my $self = shift;
    return join ',', $self->_rights;
}

sub is_trivial {
    my $self = shift;
    return $self->num_left == 1 || $self->num_right == 1;
}

sub order {
    my $self = shift;
    my $other = shift;
    return -1 if $self->num_left < $other->num_left;
    return  1 if $self->num_left > $other->num_left;
    for my $n ( 0 .. $self->num_left - 1 ) {
        return -1 if $self->left->[$n] lt $other->left->[$n];
        return  1 if $self->left->[$n] gt $other->left->[$n];
    }
    return 0;
}

sub is_compatible {
    my $self = shift;
    my $other = shift;

    my $s = $self->is_compatible_dir($other, 'l');
    my $o = $other->is_compatible_dir($self, 'l');
    
    #printf "S%1d O%1d\n", $s||0,$o||0;

    return $s||$o;
}

sub is_compatible_dir {
    my $self = shift;
    my $other = shift;
    my $dir = shift;

    $dir = $dir =~ /l(?:eft)?/i  ? '_lefts'  :
           $dir =~ /r(?:ight)?/i ? '_rights' :
           die "Unkown direction <$dir>";

    my %this_lookup = map { $_, 1 } $self->$dir;
    my $found     = 0;
    my $not_found = 0;
    for my $l ( $other->$dir ) {
        if ( exists $this_lookup{ $l } ) {
            $found++;
            return '' if $not_found; # Fail fast 1
        }
        else {
            $not_found++;
            return '' if $found; # Fail fast 2
        }
    }
    if ( $found && $not_found ) {
        # Incompatible splits
        return '';
    }
    return 1;
}

sub is_identical {
    my $self = shift;
    my $other = shift;

    # This is the fastest version so far
    return $self->stringify eq $other->stringify;

    # Fastest possible failure
    return '' if $self->num_left != $other->num_left;

    my $this_left = $self->left;
    my $that_left = $other->left;
    for my $idx ( 0 .. $self->num_left - 1 ) {
        return '' if $this_left->[$idx] ne $that_left->[$idx];
    }
    return 1;
}

__PACKAGE__->meta->make_immutable;

1;
