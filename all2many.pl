#!/usr/bin/perl
#
# Usage information updated January 21, 2005 - JDW


unless (@ARGV == 2) {
  print <<EOH;
Usage: $0  input.file  min_size_contig

    Separates a set of FASTA-format sequences in the file named
    as the first argument into individual files, each of which is
    at least as long as the  second argument. Sequences shorter
    than the second argument are ignored.  The name of each new
    file is the name of the contig.

Examples:

$0  fasta.screen.contigs  1000

creates individual files from the input file fasta.screen.contigs
for all contigs at least 1000 bases long.

$0  another.fasta.file  1

creates individual files from the input file another.fasta.file
for all contigs.

EOH
exit 0;
}

sub dump_seq {
	my($name, $seq) = @_;
	$name1 = $name;
	$name1 =~ s/>//;
        $name1 =~ s/^(\S*).*$/$1/;
	open HUNK,">$name1.fa" or die $!;
	print HUNK "$name \n";
	print HUNK $seq;
	close HUNK;
}

open BIG,$ARGV[0] or die $!;
while (<BIG>) {
  if ($_ =~ /Contig|^>/) {
 	if ($len >= $ARGV[1]) {
           dump_seq($name, $seq);
        }
     chomp;
     $name = $_;
     $seq = ''; $len = 0;
  } else {
     $seq .= $_;
     $len += length($_) - 1;
  } 
}
dump_seq($name, $seq);
close BIG;
