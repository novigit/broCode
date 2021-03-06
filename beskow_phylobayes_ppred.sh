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
generations=$2
every=$3
untils=$4
threads=$5
chain=$6

# print parameters (for debugging)
echo "Burnin: $burnin Every: $every Threads: $threads Until: $untils Chains: $chain" > beskow_phylobayes_ppred_$chain.debug.out

# load module
module load phylobayes/1.8

outname=ppred_b${burnin}_g${generations}_e${every}_u${untils}_t${threads}

## do posterior predictive test -per chain-
# do posterior predictive tests (maximum square deviation between global and taxon-specific empirical frequencies)
srun -n $threads readpb_mpi -ppred -comp -x $burnin $every $untils $chain 2> ${outname}_comp_${chain}.out
# do posterior predictive tests (mean number of distinct amino acids per site)
srun -n $threads readpb_mpi -ppred -div -x $burnin $every $untils $chain 2> ${outname}_div_${chain}.out
# compute site-specific marginal likelihoods
#srun -n $threads readpb_mpi -sitelogl -x $burnin $every $untils $chain 2> ${outname}_sitelogl_${chain}.out

# make out directory and move results there
mkdir $outname 
mv $outname* $outname/
rm *.tree
mv *.comp $outname/
mv *.div $outname/

exit

# POD Documentation
<<=cut

=head1 NAME

beskow_phylobayes_ppred.sh - Do posterior predictive checks of phylobayes chains that were run on PDC Beskow

=head1 USAGE

=over

=item beskow_phylobayes_ppred.sh <burnin> <generations> <every> <until> <threads> <chain1> <chain2> <chain3> <etc> 

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
