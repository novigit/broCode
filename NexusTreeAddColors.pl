#!/usr/bin/perl -w
use strict;
use Getopt::Long;

=head1 NAME

NexusTreeAddColors.pl - color nexus trees, add group info to labels

=head1 SYNOPSIS

NexusTreeAddColors.pl -i|--infile <sumtree.nxs> -o|--outfile <for_figtree.nxs> [-f|--figtree <figtree_options>] -p|--patterns <pattern_file> [-b|--branch-colors<file>] [-c|--color-per-group <file>] [-n|--dont-add-group-names] [-e|--exact-match]

=head1 INPUT

=head2 -i|--infile, -o|--outfile

In and out files, in nexus format

=head2 -f|--figtree <figtree_options>

A figtree option block, appended at the end of the nexus file. If not specified, a default one is provided

=head2 -p|--patterns <pattern_file>

A tabular file to color names with a color. The color is changed if the name starts with the pattern. By default, the folowing scheme is used. It is possible to use reserved words blue, green, yellow, orange, red, grey, black, purple, teal, brown,orange that will be converted to numbers:

  16777012  BH    # blue
  16777012  BQ
  16738048  m02   # green
  16738048  m07a
  16738048  BBb
  3355648   BG    # yellow
  3355648   BT
  26368     BVtw  # orange
  26368     BVwin
  3407872   BAnh1 # red

Or:

  blue   BH
  yellow BAnh1
...

Alternatively, a group name can be added in a third column, and used to rename the tags.

If the option -c|--color-per-group is activated, then the two columns must be label and group.

The group replacing option triggers exact pattern matching.

=head2 -b|--branch-colors <file>

A tab file as above. Colors terminal branches instead of labels.

=head2 -c|--color-per-group <file>

If this option is activated, groups as defined in -p|--pattern-file will be colored according to this file. Only the two first columns of the pattern file will be considered. 

Two columns: color and group:

  #8dd3c7  Acidobacteria
  #ffffb3  Alphaproteobacteria
...

=head2 -n|--dont-add-group-names

Do not add group names in front of labels.

=head2 -e|--exact-match

If set, only exact matching will be performed, instead of

  m/(\d+\s+)?$pattern(\S*)/

=head1 AUTHOR

	Katarzyna Zaremba, Katarzyna.Zaremba@icm.uu.se
        Lionel Guy, lionel.guy@icm.uu.se

=cut


my $sumtree_in;
my $figtree_out;
my $figtree_opt;
my $patterns_file;
my $branch_file;
my $color_file;
my $dont_replace;
my $exact_match;

GetOptions(
    'i|infile=s'  => \$sumtree_in,
    'o|outfile=s' => \$figtree_out,
    'f|figtree=s' => \$figtree_opt,
    'p|patterns=s'=> \$patterns_file,
    'b|branch-color=s' => \$branch_file,
    'c|color-per-group=s' => \$color_file,
    'n|dont-add-group-names' => \$dont_replace,
    'e|exact-match' => \$exact_match,
);

unless ( $sumtree_in ) {
    exec( 'perldoc', $0 );    # exec exits this program automatically
}

# Force exact match if group
$exact_match++ if $color_file;

# Constants (color codes are based on FigTree 1.4.0 "Save as.." NEXUS output)
my %def_colors = ('blue'   => 16777012,
		  'green'  => 16738048,
		  'yellow' => 2962887,
		  'orange' => 2211071,
		  'red'    => 3407872,
		  'grey'   => 8355712,
		  'black'  => 16449536,
		  'purple' => 3599223,
		  'teal'   => 16720385,
		  'brown'  => 10014976,
		  'pink'   => 65332,
		  'army'   => 10066432,
		  'sky'    => 10066177,
		  'storm'  => 16751002,
		  'salmon' => 39322,
		  'lind'   => 6684775,
		  'wood'   => 3368704,
	      );

# Parse color and groups if necessary
my %color_per_group;
if ($color_file){
  open COLORS, '<', $color_file or die "Could not open $color_file: $!\n";
  while (<COLORS>){
    chomp;
    my ($color, $group) = split;
    $color_per_group{$group} = $color;
  }
}

# Parse patterns
my %patterns;
if ($patterns_file){
  open PATT, "<", $patterns_file or die "Could not open $patterns_file\n";
  foreach (<PATT>){
    chomp;
    #next if /^#/;
    my ($color, $pattern, $group);
    my @a = split(/\s+/, $_);
    # color is only present unless color file is defined
    unless ($color_file){
      $color = shift @a;
      $color = $def_colors{$color} if $def_colors{$color};
    }
    ($pattern, $group) = @a;
    #print "$color => $pattern\n";
    if ($group) {
      $patterns{$pattern}{'group'} = $group;
      if ($color_file){
	$color = $color_per_group{$group};
	die "Color for group $group undefined\n" unless $color;
      }
      $exact_match++;
    }
    $patterns{$pattern}{'color'} = $color if $color;
    die "At least color OR group have to be defined at $_\n" 
      unless ($color or $group);
  }
}

# Parse branch colors
my %branches;
if ($branch_file){
  open BRANCH, '<', $branch_file or die "Could not open $branch_file";
  while (<BRANCH>){
    chomp;
    my ($color, $name) = split;
    $color = $def_colors{$color} if $def_colors{$color};
    $branches{$name} = $color;
  }
}

open my $IN, '<', $sumtree_in 	or die "Could not open $sumtree_in";
my $OUT;
if ($figtree_out){
  open $OUT, '>', $figtree_out or die "Could not open $figtree_out";
}
else {
  $OUT = \*STDOUT;
}

my $section = "";
while ( my $line = <$IN> ) {
  next if ( $line =~ /^\[ .+ \]$/ );
  if ($line =~ /^BEGIN (\w+);/i){
    $section = $1;
  }
  if ($line =~ /^END;/i){
    $section = "";
  }
  my $replaced;
  if ($section eq "TAXA"){
    foreach my $pattern (sort {$b cmp $a} keys %patterns){
      my $insert = "";
      $insert = "#-" unless $patterns{$pattern}{'color'} =~ /^#/;
      # If we do pattern-replacing: disable group replacing
      if (!$exact_match && $line =~ /^(\s+(\d+\s+)?$pattern\S*)/) {
	print $OUT $1 . "[&!color=" . $insert . ($patterns{$pattern}{'color'}) 
	  . "]\n";
	$replaced++;
	last;
      # else do exact match replacing and group thingy
      } elsif ($line =~ /^(\s+)$pattern/) {
	# Add group and two underscores if group defined
	print $OUT $1; 
	print $OUT $patterns{$pattern}{'group'}, "__" 
	  if ($patterns{$pattern}{'group'} && !$dont_replace);
	# Replace pattern
	print $OUT $pattern . "[&!color=" . $insert . 
	  ($patterns{$pattern}{'color'}) . "]\n";
	$replaced++;
	last;
      }
    }
    unless ($replaced){
      print $OUT "$line\n";
    }
  }
  elsif ($section eq "TREES"){
    # apply branch colors if requested: add [&!color=#-2874302]
    if ($line =~ /^\s+TREE/){
      # foreach my $name (keys %branches){
      # 	my $insert = "";
      # 	$insert = "#-" unless $branches{$name} =~ /^#/;
      # 	my $replacement = $name. '[&!color=' . $insert . $branches{$name} .
      # 	  ']';
      # 	$line =~ s/$name/$replacement/g;
      # }
      # Replacement with group names
      foreach my $org (sort {$b cmp $a} keys %patterns){
	my $group = $patterns{$org}{'group'};
	if ($group){
	  my $replacement = $group . "__" . $org;
	  $line =~ s/$org/$replacement/ 
	    if ($patterns{$org}{'group'} && !$dont_replace);
	  # print $OUT $patterns{$org}{'group'}, "__" 
	  #   if ($patterns{$org}{'group'} && !$dont_replace);
	}
      }
    }
    if ($line =~ /(TREE .* = (\[&U\] )?)(.+)/) {
      my $desc = $1;
      my $newick = $3;
      $newick =~ s/\)([\d.]+):/)[&bootstrap=$1]:/g;
      $newick =~ s/\)([\d.]+);$/)[&bootstrap=$1];/;
      print $OUT "$desc $newick\n";
    }
    # specific for trees that earlier looked at and changed
    elsif ($line =~ /tree 0 = \[&U\] (.+)/) {
      my $newick = $1;
      $newick =~ s/\)\[&label=([\d.]+)\]:/)[&bootstrap=$1]:/g;
      $newick =~ s/\)\[&label=([\d.]+)\];$/)[&bootstrap=$1];/;
      print $OUT "TREE 0 = \[&U\] $newick\n";
    }
    # change bootstrap
    elsif ($line =~ /(^.*tree Bioperl_.* = \[&R\]) (.+)/) {
      my $desc = $1;
      my $nexus = $2;
      $nexus =~ s/\[(0?\.?[\d.]+)\]/[&bootstrap=$1]/g;
      #$nexus =~ s/\[\]/)/g;
      print $OUT "$desc $nexus\n";
    }
    else {
      print $OUT "$line\n";
    }
  }
  # remove previous figtree formatting
  elsif ($section eq "figtree") {
    # do nothing
  }
  else {
    print $OUT $line unless $replaced;
  }
}

close($IN);
if ($figtree_opt){
  close($OUT);
  system("cat $figtree_opt >> $figtree_out");
}
else {
  my $str = figTreeOpts();
  print $OUT $str;
  close($OUT);
}


sub figTreeOpts {
    my $str = <<'EOF';
begin figtree;
        set appearance.backgroundColour=#-1;
        set appearance.branchColorAttribute="User Selection";
        set appearance.branchLineWidth=1.0;
        set appearance.foregroundColour=#-16777216;
        set appearance.selectionColour=#-2144520576;
        set branchLabels.displayAttribute="bootstrap";
        set branchLabels.fontName="Arial";
        set branchLabels.fontSize=10;
        set branchLabels.fontStyle=0;
        set branchLabels.isShown=true;
        set branchLabels.significantDigits=4;
        set layout.layoutType="RECTILINEAR";
        set layout.zoom=0;
        set nodeBars.barWidth=4.0;
        set nodeLabels.displayAttribute="Node ages";
        set nodeLabels.fontName="Arial";
        set nodeLabels.fontSize=8;
        set nodeLabels.fontStyle=0;
        set nodeLabels.isShown=false;
        set nodeLabels.significantDigits=4;
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
        set scaleAxis.majorTicks=0.1;
        set scaleAxis.origin=0.0;
        set scaleAxis.significantDigits=4;
        set scaleBar.automaticScale=true;
        set scaleBar.fontSize=10.0;
        set scaleBar.isShown=true;
        set scaleBar.lineWidth=1.0;
        set scaleBar.scaleRange=0.0;
        set scaleBar.significantDigits=4;
        set tipLabels.displayAttribute="Names";
        set tipLabels.fontName="Arial";
        set tipLabels.fontSize=12;
        set tipLabels.fontStyle=0;
        set tipLabels.isShown=true;
        set tipLabels.significantDigits=4;
        set trees.order=false;
        set trees.orderType="increasing";
        set trees.rootingType="User Selection";
        set trees.transform=false;
        set trees.transformType="cladogram";
end;
EOF
    return $str;
}
