#!/bin/bash -l

script=$0
# state usage in POD format
function usage() {
    pod2usage $script >&2
    exit -1
}

# check if arguments are provided
if (( $# < 4 )); then
    echo "Not enough arguments" >&2;
    usage;
fi

# set options
phylip=$1
model=$2
chain=$3
mode=$4

# load tools
module load phylobayes/1.8

# start phylobayes chain
# by NOT setting -S, and saving every generation, we allow for posterior predictive tests
# setting -dc removes constant columns from the alignment and improves mixing (see manual)
if [ "$mode" = "start" ]; then
    srun -n 32 pb_mpi -d $phylip -cat -$model $chain &
fi

# continue phylobayes chain
if [ "$mode" = "continue" ]; then
    srun -n 32 pb_mpi $chain &
fi

# stop chains before time runs out
sleep 23h && sleep 50m
echo 0 > $chain.run

exit

# POD Documentation
<<=cut

=head1 NAME

beskow_phylobayes.sh - Run one phylobayes chain on Beskow

=head1 USAGE

beskow_phylobayes.sh <phylip> <model> <chain> <mode>

=head1 DESCRIPTION

This script is called by the script beskow_submit_phylobayes.sh. 

It tells Beskow to run one phylobayes chain on the alignment <phylip> with the model <model>, for <hours> hours

It has two modes: start and continue. Start initates the chain, while continue continues a chain.

Currently hard-coded to running 128 cores per chain, for 23 hours and 50 minutes

=head1 NOTE

This script should be used only on Beskow, because its large CPU power.

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)

=cut
