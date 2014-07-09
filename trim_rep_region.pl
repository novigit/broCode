#!/usr/bin/perl

# scripting by Guillaume Reboul : reboul_guillaume@yahoo.fr

###############
#This script trim the repeat region found by repeat-match and create a new fasta file for the trimmed sequence
#Plus, it print some general information about the repeat region which can be put in q file with >> <output info file>
#Can be used with shell loop but not with multifasta file
###############

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;

#quit and print the usage if there is no 3 arguments
die "usage : trim_rep_region.pl <repeat-match file> <fasta in> <fasta out>\n" if @ARGV!=3;



###############
# global variables
###############

my ($file1,$file2,$file3)=@ARGV;

open(RepFile,$file1);
my $fastaIn=Bio::SeqIO->new(-file=>$file2, -format=>"fasta");
my $fastaOut=Bio::SeqIO->new(-file=>">".$file3, -format=>"fasta");

my $seqFastaIn = $fastaIn->next_seq;
my $idFastaIn=$seqFastaIn->id;
my $sequenceFastaIn=$seqFastaIn->seq;

my $check=0;
my %repeat;



###############
# repeat storage in the hash %repeat
###############

while(<RepFile>){
    if ($_ =~ /\s+(\d+)\s+(\d+)r\s+(\d+)/){
        $check++;
	$repeat{$check}{'startF'}=$1;
	$repeat{$check}{'startR'}=$2;
	$repeat{$check}{'repLen'}=$3;
    }
}
close(RepFile);



###############
# main treatment
###############

#if there is no repeat region
if ($check == 0){
    $fastaOut->write_seq($seqFastaIn);
    die "No reverse repeat-matches in the repeat-match file $file1 associate to the fasta file $file2\n";
}

#if there is one repeat region
if ($check == 1){
    #if the repeat region begin at the first nucleotide (+-5)
    if ($repeat{'1'}{'startF'} <= 1+5){
	trim_seq($sequenceFastaIn,$repeat{'1'}{'repLen'},'f');
    }
    #if the repeat region gegin at the last nucleotide (+-5)
    elsif ($repeat{'1'}{'startR'} >= length($sequenceFastaIn)-5){
	trim_seq($sequenceFastaIn,$repeat{'1'}{'repLen'},'r');
    }
    #else, the repeat region is in the sequence so untreated here
    else {
	$fastaOut->write_seq($seqFastaIn);
	#print message if the repeat region is more than 50 bases
	print "untreated : $idFastaIn : $repeat{'1'}{'repLen'}\n" if ($repeat{'1'}{'repLen'} >= 50);
        die "No reverse repeat-matches in the repeat-match file $file1 associate to the fasta file $file2\n";
    }
    print "treated : $idFastaIn : $repeat{'1'}{'repLen'}\n";
}
#if there are more than 1 repeat region
else{
    my @finalGoodOrder;
    my $finalLengthToTrim=0;
    my $strandGoodOrder='';
    my $currentLengthToTrim;
    #for all the repeat in the %repeat
    for (1..$check){
	#if the repeat region begin at the first nucleotide (+-5)
	if ( $repeat{$_}{'startF'} <= 1+5 ){
	    #the key of the repeat is stocked
	    my @testCurrentRow=$_;
	    #call the recursive fonction
	    @testCurrentRow=order_row($_,$check,'f',\@testCurrentRow);
	    #calculation of the entire repeat length
	    $currentLengthToTrim=0;
	    foreach (@testCurrentRow){
		$currentLengthToTrim+=$repeat{$_}{'repLen'};
	    }
	    #test to keep the longest repeat found
	    if ($currentLengthToTrim >= $finalLengthToTrim){
		#if the latest repeat was more than 50 base, the user is noticed
		print "untreated : $idFastaIn : $finalLengthToTrim\n" if ($finalLengthToTrim >= 50);
                $strandGoodOrder='f';
                $finalLengthToTrim=$currentLengthToTrim;
		@finalGoodOrder=@testCurrentRow;
	    }
	    #if the newest repeat is more than 50 base but less than the stocked repeat, the user is noticed too
	    else {
                print "untreated : $idFastaIn : $currentLengthToTrim\n" if ($currentLengthToTrim >= 50);
            }
	}
	#if the repeat region gegin at the last nucleotide (+-5)
	elsif ($repeat{$_}{'startR'} >= length($sequenceFastaIn)-5){
	    my @testCurrentRow=$_;
            @testCurrentRow=order_row($_,$check,'r',\@testCurrentRow);
	    $currentLengthToTrim=0;
	    foreach (@testCurrentRow){
                $currentLengthToTrim+=$repeat{$_}{'repLen'};
            }
            if ($currentLengthToTrim >= $finalLengthToTrim){
		print "untreated : $idFastaIn : $finalLengthToTrim\n" if ($finalLengthToTrim >= 50);
                $strandGoodOrder='r';
                $finalLengthToTrim=$currentLengthToTrim;
                @finalGoodOrder=@testCurrentRow;
            }
	    else {
		print "untreated : $idFastaIn : $currentLengthToTrim\n" if ($currentLengthToTrim >= 50);
	    }
	}
    }
    # to treat only the sequence which begin or finish by a repeat region
    if ($finalLengthToTrim != 0 && $strandGoodOrder ne ''){
	print "treated : $idFastaIn : $finalLengthToTrim\n";
	trim_seq($sequenceFastaIn,$finalLengthToTrim,$strandGoodOrder);
    }
    else {
	$fastaOut->write_seq($seqFastaIn);
	print "untreated : $idFastaIn : $repeat{'1'}{'repLen'}\n" if ($repeat{'1'}{'repLen'} >= 50);
	die "No reverse repeat-matches in the repeat-match file $file1 associate to the fasta file $file2\n";
    }
}



###############
# fonctions
###############

# the recursive fonction to find the order making up the entire repeat region
sub order_row {
    my ($previousGoodRow,$nbRow,$strand,$sortRow)=@_;
    if ($strand eq 'f'){
	#for all the repeat in the %repeat
	for (my $testRow=$nbRow;$testRow>0;$testRow--){
	    #test if the tested repeat is the follower of the repeat previously stocked
	    if ($repeat{$testRow}{'startF'}-$repeat{$previousGoodRow}{'repLen'}-$repeat{$previousGoodRow}{'startF'} == 1){
		#if it is, it became the new stocked repeat
		push(@{$sortRow},$testRow);
		#and the fonction is called with it to found the next repeat region
		@{$sortRow}=order_row($testRow,$nbRow,'f',\@{$sortRow});
	    }
	}
    }
    elsif ($strand eq 'r'){
	for (my $testRow=$nbRow;$testRow>0;$testRow--){
            if ($repeat{$previousGoodRow}{'startR'}-$repeat{$previousGoodRow}{'repLen'}-$repeat{$testRow}{'startR'} == 1){
                push(@{$sortRow},$testRow);
		@{$sortRow}=order_row($testRow,$nbRow,'r',\@{$sortRow});
            }
        }
    }
    #return the table with the sorted number of the key of the %repeat making up the entire repeat region
    return @{$sortRow};
}

#fonction to trim the fasta sequence and print the new one in the new file defined by the user
sub trim_seq {
    my ($seqToTrim,$lengthToTrim,$strand)=@_;
    my $newSeq;
    if ($strand eq 'f'){
	$newSeq=Bio::Seq->new(-id=>$idFastaIn, -seq=>substr($seqToTrim,$lengthToTrim));
    }
    elsif ($strand eq 'r'){
        $newSeq=Bio::Seq->new(-id=>$idFastaIn, -seq=>substr($seqToTrim,0,length($seqToTrim)-$lengthToTrim));
    }
    else {
	die "problem\n";
    }
    $fastaOut->write_seq($newSeq);
}
