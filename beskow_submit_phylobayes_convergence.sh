#!/bin/bash

script=$0
# state usage
function usage() {
    pod2usage $script >&2
    exit -1
}
function help_doc() {
    pod2text $script >&2
    exit -1
}

# check number of arguments
if (( $# == 0 || $# >= 2 && $# <= 4 )); then
    echo -e "Error: not enough arguments" "\n" >&2;
    usage;
fi

# check if invoked on Tegner
if [[ $(hostname -s) != tegner* ]]; then
    echo -e "Error: $script can only be run on Tegner" "\n" >&2;
    help_doc;
fi

# state options
while getopts ":b:o:h" opt; do
    case $opt in
	b) burnin=${OPTARG};;
	o) outname=${OPTARG};;
	h) help_doc;;
	*) usage;;
    esac
done

# remove stated options from $*
shift $(( OPTIND - 1 )) # OPTIND's value is the index of the next argument to be processed
chains=$*

# submit job
sbatch -A 2020-5-473 -J $outname -o $outname.o -e $outname.e -t 20:00 --nodes=1 --ntasks-per-node=1 --mail-type=BEGIN,END,FAIL beskow_phylobayes_convergence.sh $burnin $outname $chains

# POD Documentation
<<=cut

=head1 NAME

beskow_submit_phylobayes_convergence.sh - Submit job to check convergence of phylobayes chains that were run on PDC Beskow

=head1 USAGE

beskow_submit_phylobayes_convergence.sh -b <burnin> -o <outname> <chain1> <chain2> <etc>

<chains> have to be last arguments!

<outname> will be the name of the output directory and the prefix for all outfiles

For more documentation:

beskow_submit_phylobayes_convergence.sh -h

=head1 DESCRIPTION

Calls the script beskow_phylobayes_convergence.sh to submit convergence analyses of phylobayes chains.

=head1 NOTE

This script should be used only on Tegner, because Tegners purpose is for post-processing.

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)

=cut
