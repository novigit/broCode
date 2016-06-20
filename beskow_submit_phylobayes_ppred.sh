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
# if [[ $(hostname -s) != tegner* ]]; then
#     echo -e "Error: $script can only be run on Tegner" "\n" >&2;
#     help_doc;
# fi

# state options
while getopts ":b:e:t:h" opt; do
    case $opt in
	b) burnin=${OPTARG};;
	e) every=${OPTARG};;
	t) threads=${OPTARG};;
	h) help_doc;;
	*) usage;;
    esac
done

# remove stated options from $*
shift $(( OPTIND - 1 )) # OPTIND's value is the index of the next argument to be processed
chains=$*

# print parameters (for debugging)
echo "Burnin: $burnin Every: $every Threads: $threads Chains: $chains" > beskow_submit_phylobayes_ppred.debug.out

outname=ppred_b${burnin}_e${every}

# submit job
sbatch -A m.2015-1-273 -J $outname -o $outname.o -e $outname.e -t 8:00:00 --nodes=1 --ntasks-per-node=$threads --mail-type=BEGIN,END,FAIL beskow_phylobayes_ppred.sh $burnin $every $threads $chains

# POD Documentation
<<=cut

=head1 NAME

beskow_submit_phylobayes_ppred.sh - Submit job to perform posterior predictive checks of phylobayes chains that were run on PDC Beskow

=head1 USAGE

beskow_submit_phylobayes_ppred.sh -b <burnin> -e <every> -t <threads> <chain1> <chain2> <etc>

<chains> have to be last arguments!

<outname> will be the name of the output directory and the prefix for all outfiles

For more documentation:

beskow_submit_phylobayes_ppred.sh -h

=head1 DESCRIPTION

Calls the script beskow_phylobayes_ppred.sh to submit posterior predictive checks of phylobayes chains.

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)

=cut
