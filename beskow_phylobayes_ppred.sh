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

# check if arguments are provided
if (( $# < 3 )); then
    echo "Not enough arguments" >&2;
    usage;
fi

# options
burnin=$1
shift
every=$1
shift
threads=$1
shift
chains=$*

# print parameters (for debugging)
echo "Burnin: $burnin Every: $every Threads: $threads Chains: $chains" > beskow_phylobayes_ppred.debug.out

# load module
module load phylobayes/1.5a

outname=ppred_b${burnin}_e${every}_t${threads}
# do posterior predictive test -per chain-
for chain in $chains; do

    ## do posterior predictive tests (mean number of distinct amino acids per site)
    # aprun -n $threads readpb_mpi -ppred -div -x $burnin $every $chain 2> ${outname}_div_${chain}.out
    # do posterior predictive tests (maximum square deviation between global and taxon-specific empirical frequencies)
    aprun -n $threads readpb_mpi -ppred -comp -x $burnin $every $chain 2> ${outname}_comp_${chain}.out

done

# make out directory and move results there
mkdir $outname
mv $outname* $outname

exit

# POD Documentation
<<=cut

=head1 NAME

beskow_phylobayes_ppred.sh - Do posterior predictive checks of phylobayes chains that were run on PDC Beskow

=head1 USAGE

=over

=item beskow_phylobayes_ppred.sh <burnin> <threads> <chain1> <chain2> <chain3> <etc> 

=item <chains> have to be last arguments!

=back

=head1 DESCRIPTION

This script will simulate datasets based on parameters sampled from the posterior probability distribution and calculate per simulated dataset a test statistic. In this case the test statistics are the mean number of distinct amino-acids per site, and the maximum square deviation between global and taxon-specific empirical frequencies

=head1 NOTE

=over

=item It seems that readpb_mpi does not work properly on Tegner. Try on beskow instead

=item This script is called by the script beskow_submit_phylobayes_ppred.sh and is the only way this script works

=back

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)

=cut