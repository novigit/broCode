#!/usr/bin/env perl
package remove_species_in_splitSet;
use strict;
use warnings;
use feature ':5.10';
use Moose;
use MooseX::Types::Path::Class;
with 'MooseX::Getopt';

use FindBin;
use lib "$FindBin::Bin";
use Phylogeny::Tree::SplitSet;

has 'file' => (
    is => 'ro',
    isa => 'Path::Class::File',
    coerce => 1,
    required => 1,
    documentation => 'File of splits',
);

has outfile => (
    is => 'ro',
    isa => 'Path::Class::File',
    coerce => 1,
    required => '',
    documentation => 'Outfile to put result in',
);

has species => (
    is            => 'ro',
    isa           => 'Str',
    required      => '',
    documentation => 'Species to remove from the split set, separated by commas',
);


sub run {
    my $self = shift;

    my $split_set = $self->get_split_set( $self->file );
    my $species = $self->get_species( $self->species );

    my $new_split_set = $split_set->remove_species(@{ $species });

    my $str = $new_split_set->stringify;
    my $OUT = $self->outfile->openw;
    say $OUT $str;
}

sub get_split_set {
    my $self = shift;
    my $splits_file = shift;

    open my $IO_TREE, '<', $splits_file or die "Can't open file $splits_file: $!";
    my $str = join '', <$IO_TREE>;

    my $split;
    eval {
        $split = Phylogeny::Tree::SplitSet->from_string( $str );
    }; if ( $@ ) {
        die "Can't parse split from $splits_file: $@";
    }
    return $split;
}

sub get_species {
    my $self = shift;
    my $species = shift;

    my @species = split(/,/, $species);
    return \@species;
}


remove_species_in_splitSet->new_with_options->run unless caller;

__PACKAGE__->meta->make_immutable;

1;
