#!/bin/bash -l

script=$0
# state usage in POD format
function usage() {
    pod2usage $script >&2
    exit -1
}

# state help
function help_doc() {
    pod2text $script >&2
    exit -1
}

# check number of arguments
if (( $# == 0 || $# >= 2  && $# <= 6 )); then
    echo -e "Error: not enough arguments" "\n";
    usage;
fi

# state options
while getopts ":p:m:s:h" opt; do
    case $opt in
	p) phylip=${OPTARG};;
	m) model=${OPTARG};;
	s) setting=${OPTARG};;
	h) help_doc;;
	*) usage;;
    esac
done
shift $(( OPTIND - 1 ))
chains=$*

# submit job
jobname=$(basename $phylip .phylip)
for chain in $chains;
do
    # start chains
    if [ "$setting" = "start" ]; then
	mkdir chain_backup
	sbatch -A 2020-5-473 -J $jobname -o $jobname.o -e $jobname.e -t 23:59:00 --nodes=1 --ntasks-per-node=32 --mail-type=BEGIN,END,FAIL beskow_phylobayes.sh $phylip $model $chain start
    fi

    # continue chains
    if [ "$setting" = "continue" ]; then
	rsync -av --partial --inplace --append --progress ${chain}*.* chain_backup/
	sbatch -A 2020-5-473 -J $jobname -o $jobname.o -e $jobname.e -t 23:59:00 --nodes=1 --ntasks-per-node=32 --mail-type=BEGIN,END,FAIL beskow_phylobayes.sh $phylip $model $chain continue
    fi
done

exit

# POD Documentation
<<=cut

=head1 NAME

beskow_submit_phylobayes.sh - Submit job to run phylobayes chains on PDC Beskow

=head1 USAGE

beskow_submit_phylobayes.sh -p <phylip> -m <model> -s <setting> <chain1> <chain2> <etc>

<chains> have to be last arguments!

For more documentation:

beskow_submit_phylobayes.sh -h

=head1 DESCRIPTION

Calls the script beskow_phylobayes.sh to submit phylobayes chains on PDC Beskow

It has two modes: start and continue. Start to initiate phylobayes chains and continue to continue chains.

Currently hard-coded to running 128 cores per chain, for 23 hours and 50 minutes

=head1 EXAMPLES

beskow_submit_phylobayes.sh -p alignment.phylip -m gtr -s start chain1 chain2 chain3 chain4

beskow_submit_phylobayes.sh -p alignment.phylip -m gtr -s continue chain1 chain2 chain3 chain4

=head1 NOTE

This script should be used only on PDC Beskow, because of its large CPU power.

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)

=cut
