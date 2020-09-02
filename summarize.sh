#!/bin/bash -l
#SBATCH -A m.2015-1-273
#SBATCH -J summarize
#SBATCH -t 30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -e summarize_error.e
#SBATCH -o summarize_output.o

burnin=10000
subs=10
base=sr4

files=

module load phylobayes/1.5a

lines=$(wc -l ${base}_1.trace)
gen=$(echo $lines | cut -f1 -d' ')
summary_folder=${base}_$gen

aprun -n 1 tracecomp -x $burnin $files
summary=${base}_$gen.trace_$burnin.summary
cat summarize_error.e > $summary
cat summarize_output.o >> $summary

aprun -n 1 bpcomp -x $burnin $subs $files
sed -e "s/)1:/)1.0:/g" bpcomp.con.tre > bpcomp.con.tree

for i in 1 2 3 4
do
        aprun -n 1 bpcomp -x $burnin $subs ${base}_$i
        sed -e "s/)1:/)1.0:/g"  ${base}_$i.con.tre >  ${base}_$i.con.tree
done

mkdir $summary_folder
cp bpcomp.* $summary_folder
cp $base*tree $summary_folder



