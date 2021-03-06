#!/bin/bash -l

script=$0
# state usage in POD format
function usage() {
    pod2usage $script >&2
    exit -1
}

# check if arguments are provided
if (( $# < 7 )); then
    echo "Not enough arguments" >&2;
    usage;
fi

# set options
fwdReads=$1
rvrReads=$2
unpReads=$3
threads=$4
outpdir=$5
memlimit=$6
kmers=$7

# load tools
module load SPAdes/3.8.1

# start spades assembly
spades.py --sc -1 $fwdReads -2 $rvrReads -s $unpReads -t $threads -o $outpdir -m $memlimit -k $kmers

exit

# POD Documentation
<<=cut

=head1 NAME

tegner_spades.sh - Run spades assembly on tegner

=head1 USAGE

tegner_spades.sh <fwdReads> <rvrReads> <unpReads> <threads> <outpdir> <memlimit> <kmers>

=head1 DESCRIPTION

Run spades assembly on tegner on pdc

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)

=cut
