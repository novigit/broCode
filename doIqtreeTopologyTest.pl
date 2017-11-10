#!/usr/bin/perl

# libraries
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

# options
my ($constraints, $bootstraptrees, $alignment, $model, $threads, $guidetree); #input
my ($results, $outputdir); # output
my $help;
GetOptions(
    'c|constraints=s'        => \$constraints,
    'b|bootstraptrees=s'     => \$bootstraptrees,
    'a|alignment=s'          => \$alignment,
    'm|model=s'              => \$model,
    't|threads=i'            => \$threads,
    'g|guidetree=s'          => \$guidetree,
    'r|results=s'            => \$results,
    'h|help'                 => \$help,
    );

pod2usage(-exitval => 1, -verbose => 2) if $help;
pod2usage(-msg => "ERROR: Please provide a constraints file" ) unless $constraints;
pod2usage(-msg => "ERROR: Please provide bootstrap trees" ) unless $bootstraptrees;
pod2usage(-msg => "ERROR: Please provide an alignment" ) unless $alignment;
pod2usage(-msg => "ERROR: Please provide a substitution model" ) unless $model;
pod2usage(-msg => "ERROR: Please specify number of threads" ) unless $threads;
pod2usage(-msg => "ERROR: Please specify an output results file" ) unless $results;

my %hash;

print "Guidetree detected, invoke PMSF model under ML tree searches and topology test", "\n" if $guidetree;

# CONSTRAINT ML SEARCH
# rm mlConstraints.treelist file from previous runs of this script
system("rm mlConstraints.treelist") if (-e 'mlConstraints.treelist');
# read constraint file, split into separate constraints, each with a name #
open CONSTRAINTS, "<", $constraints;
while (my $line = <CONSTRAINTS>) {
    my ($constraintName, $constraintTree) = split "\t", $line;

    my $iqtreeTreeName = "Tree $.";
    print $iqtreeTreeName, "\n";
    $hash{$iqtreeTreeName} = $constraintName;

    # create constraint file
    my $constraintFile = $constraintName . '.constraint';
    open CONSTRAINTFILE, ">", $constraintFile;
    print CONSTRAINTFILE $constraintTree;
    close CONSTRAINTFILE;

    # run constraint ML search for each specified constraint.. #
    
    # if guidetree is given:
    if ($guidetree) {
	my $cmd1_pmsf = "iqtree-omp -s $alignment -m $model -ft $guidetree -g $constraintFile -nt $threads -pre $constraintName.PMSF -quiet -redo";
	print STDERR "Searching for ML tree under PMSF model under constraint $constraintName, using $threads threads" , "\n";
	system($cmd1_pmsf);
    } else {
    # if no guidetree is given:
	my $cmd1 = "iqtree-omp -s $alignment -m $model -g $constraintFile -nt $threads -pre $constraintName -quiet -redo";
	print STDERR "Searching for ML tree under constraint $constraintName, using $threads threads" , "\n";
	system($cmd1);
    }
    # collect constraint ML treefiles
    if ($guidetree) {
	system("cat $constraintName.PMSF.treefile >> mlConstraints.treelist");
    } else {
	system("cat $constraintName.treefile >> mlConstraints.treelist");
    }
}
close CONSTRAINTS;

# POOL TREES
# pool constraint ML trees with bootstrap trees
my $cmd2 = "cat mlConstraints.treelist $bootstraptrees > forTopologyTest.treelist";
print STDERR "Pooling obtained ML trees with $bootstraptrees", "\n";
system($cmd2);

# EXECUTE TOPOLOGY TEST
my $cmd3;
if ($guidetree) {
    $cmd3 = "iqtree-omp -s $alignment -m $model -ft $guidetree -nt $threads -pre TopologyTest -z forTopologyTest.treelist -n 1 -zb 10000 -zw -au -fixbr -wsl -quiet -redo";
} else {
    $cmd3 = "iqtree-omp -s $alignment -m $model -nt $threads -pre TopologyTest -z forTopologyTest.treelist -n 1 -zb 10000 -zw -au -fixbr -wsl -quiet -redo";
}

print STDERR "Executing topology test, it may take a while", "\n";
system($cmd3);

# ORGANIZE RESULTS
open IQTREE, "<", 'TopologyTest.iqtree';
open RESULTS, ">", $results;

foreach my $tree (keys %hash) {
    print RESULTS $tree, " refers to constraint ", $hash{$tree}, "\n";
}

print RESULTS "Other trees are bootstrap trees", "\n";
my $printSwitch = 0;
while (my $line = <IQTREE>) {
    print RESULTS $line if ($printSwitch == 1);
    $printSwitch = 1 if ($line =~ m/USER TREES/);
    $printSwitch = 0 if ($line =~ m/TIME STAMP/);
}
close RESULTS;
close IQTREE;

=head1 NAME

doIqtreeTopologyTest.pl - wrapper for performing topology tests with IQ-TREE

=head1 SYNOPSIS

doIqtreeTopologyTest.pl -c|--constraints <FILE> -b|--bootstraptrees <FILE> -a|--alignment <FILE> -m|--model <SUBSTITUTION_MODEL> -t|--threads <INTEGER> [ -g|--guidetree <GUIDETREE> ] -r|--results <FILE> -o|--outputdir <OUTDIR> -h|--help

=head1 INPUT

=over

=item B<-c>, B<--constraints> I<FILE>

List of newick trees defining the constraints you want to consider

NOTE: Within the constraint, each group within two brackets must have at least 2 taxa, i.e. (taxa1,taxa2). (taxa1) will yield an error

=item B<-b>, B<--bootstraptrees> I<FILE>

List of bootstrap trees that are used to improve the topology test accuracy
Make sure that they are bootstrap trees with branch length information

=item B<-a>, B<--alignment> I<FILE>

Alignment. Should be the same that was used to generate the bootstrap trees

=item B<-m>, B<--model> I<SUBSTITUTION_MODEL>

Substitution model of evolution. Should be the same that was used to generate the bootstrap trees

=item B<-t>, B<--threads> I<INTEGER>

Number of threads

=item B<-g>, B<--guidetree> I<FILE>

Optional. Guidetree that is necessary for the PMSF approximation model

=back

=head1 OUTPUT

=over

=item B<-r>, B<--results> I<FILE>

Main results of the topology test. Excerpt from the .iqtree output file. Shows p-value table across all tests
Above the table a legend describing what is Tree1, Tree2, Tree3, etc.

=back

=head1 DESCRIPTION

Will execute a topology test for two or more considered topologies using IQ-TREE.

First executes a constraint tree searches for all provided constraints.

Then pooles obtained constraint ML trees with the bootstrap trees.

Then executes the topology test.

Finally summarizes the p-values in a table.

=head1 DEPENDENCIES

iqtree-omp

=head1 AUTHOR

Joran Martijn (L<joran.martijn@icm.uu.se>)

=head1 DATE

October 2017

=head1 COPYRIGHT AND LICENSE

Copyright 2017- Joran Martijn

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
