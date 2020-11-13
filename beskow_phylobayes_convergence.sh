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

# check dependencies
type -P plot_pb_stats.R >/dev/null 2>&1 || { echo -e "Error: plot_pb_stats.R not found" "\n" >&2 ; help_doc; }
type -P upp_convertNewick2Nexus.R >/dev/null 2>&1 || { echo -e "Error: upp_convertNewick2Nexus.R not found" "\n" >&2 ; help_doc; }

# options
burnin=$1
shift
outname=$1
shift
chains=$*

# load module
module load phylobayes/1.8 R/3.2.0

# do plot_pb_stats.R (do one plot with burnin = 0, to check if burnin value is ok)
plot_pb_stats.R folder=. burnin=0 output=${outname}.b0.Rplot.png $chains
plot_pb_stats.R folder=. burnin=$burnin output=${outname}.Rplot.png $chains

# do tracecomp
tracecomp -o $outname -x $burnin $chains
    
# calculate $every, such that you sample 200 trees
# chain1=$1
# generations=$(wc -l $chain1.trace | cut -f 1 -d ' ')
# every=$(echo "($generations - $burnin) / 200" | bc)

# do bpcomp for all chains
bpcomp -c 0.1 -o $outname -x $burnin 10 $chains
# make tree figtree-readable
sed -i -e "s/)1:/)1.0:/g" $outname.con.tre

# do bpcomp per chain
for chain in $chains; 
do 
    bpcomp -c 0.1 -o $outname.$chain -x $burnin 10 $chain; 
    sed -i -e "s/)1:/)1.0:/g" $outname.$chain.con.tre
    tracecomp -o $outname.$chain -x $burnin $chain;
done

# order tree with figtree block
# figtreeblock=$(sed '0,/^__FIGTREEBLOCK__$/d' "$0" ) # load figtreeblock from within script (see below)
# for tree in *.con.tre;
# do
#     upp_convertNewick2Nexus.R $tree
#     cat $tree.nex > $tree.nex.block
#     printf '%s\n' "$figtreeblock" >> $tree.nex.block
# done

# move outfiles to outdir
mkdir $outname
mv $outname.* $outname/

exit

# POD Documentation
<<=cut

=head1 NAME

beskow_phylobayes_convergence.sh - Check convergence of phylobayes chains that were run on PDC Beskow

=head1 USAGE

=over

=item beskow_phylobayes_convergence.sh <burnin> <outname> <chain1> <chain2> <chain3> <etc> 

=item <chains> have to be last arguments!

=item <outname> will be the name of the output directory and the prefix for all outfiles

=back

=head1 DESCRIPTION

This script will do the following:

=over

=item Run plot_pb_stats.R with burnin 0 and burnin <burnin> 

=item Run tracecomp

=item Run bpcomp on the collection of chains, and per chain

=item Make each tree ordered and viewable with FigTree

=back

=head1 DEPENDENCIES

All scripts should be within the PATH

=over

=item R and the R modules RColorBrewer and Signal

=item plot_pb_stats.R

=item upp_convertNewick2Nexus.R 

=back

=head1 NOTE

=over

=item This script should be used only on Tegner, because Tegners purpose is for post-processing.

=item This script is called by the script beskow_submit_phylobayes_convergence.sh and is the only way this script works

=back

=head1 AUTHOR

Joran Martijn (joran.martijn@icm.uu.se)

=cut

__FIGTREEBLOCK__
begin figtree;
        set appearance.backgroundColorAttribute="Default";
        set appearance.backgroundColour=#-1;
        set appearance.branchColorAttribute="User selection";
        set appearance.branchLineWidth=1.0;
        set appearance.branchMinLineWidth=0.0;
        set appearance.branchWidthAttribute="Fixed";
        set appearance.foregroundColour=#-16777216;
        set appearance.selectionColour=#-2144520576;
        set branchLabels.colorAttribute="User selection";
        set branchLabels.displayAttribute="label";
        set branchLabels.fontName="Arial";
        set branchLabels.fontSize=8;
        set branchLabels.fontStyle=0;
        set branchLabels.isShown=true;
        set branchLabels.significantDigits=4;
        set layout.expansion=411;
        set layout.layoutType="RECTILINEAR";
        set layout.zoom=0;
        set legend.attribute="label";
        set legend.fontSize=10.0;
        set legend.isShown=false;
        set legend.significantDigits=4;
        set nodeBars.barWidth=4.0;
        set nodeBars.displayAttribute=null;
        set nodeBars.isShown=false;
        set nodeLabels.colorAttribute="User selection";
        set nodeLabels.displayAttribute="Node ages";
        set nodeLabels.fontName="Arial";
        set nodeLabels.fontSize=8;
        set nodeLabels.fontStyle=0;
        set nodeLabels.isShown=false;
        set nodeLabels.significantDigits=4;
        set nodeShape.colourAttribute="User selection";
        set nodeShape.isShown=false;
        set nodeShape.minSize=10.0;
        set nodeShape.scaleType=Width;
        set nodeShape.shapeType=Circle;
        set nodeShape.size=4.0;
        set nodeShape.sizeAttribute="Fixed";
        set polarLayout.alignTipLabels=false;
        set polarLayout.angularRange=0;
        set polarLayout.rootAngle=0;
        set polarLayout.rootLength=100;
        set polarLayout.showRoot=true;
        set radialLayout.spread=0.0;
        set rectilinearLayout.alignTipLabels=false;
        set rectilinearLayout.curvature=0;
        set rectilinearLayout.rootLength=100;
        set scale.offsetAge=0.0;
        set scale.rootAge=1.0;
        set scale.scaleFactor=1.0;
        set scale.scaleRoot=false;
        set scaleAxis.automaticScale=true;
        set scaleAxis.fontSize=8.0;
        set scaleAxis.isShown=false;
        set scaleAxis.lineWidth=1.0;
        set scaleAxis.majorTicks=1.0;
        set scaleAxis.origin=0.0;
        set scaleAxis.reverseAxis=false;
        set scaleAxis.showGrid=true;
        set scaleBar.automaticScale=true;
        set scaleBar.fontSize=8.0;
        set scaleBar.isShown=true;
        set scaleBar.lineWidth=1.5;
        set scaleBar.scaleRange=0.0;
        set tipLabels.colorAttribute="User selection";
        set tipLabels.displayAttribute="Names";
        set tipLabels.fontName="sansserif";
        set tipLabels.fontSize=10;
        set tipLabels.fontStyle=0;
        set tipLabels.isShown=true;
        set tipLabels.significantDigits=4;
        set trees.order=true;
        set trees.orderType="increasing";
        set trees.rooting=true;
        set trees.rootingType="User Selection";
        set trees.transform=false;
        set trees.transformType="cladogram";
end;
