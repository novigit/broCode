#!/usr/bin/perl
#
# Version 1.0 (June 13, 2012)
# Author: Laura Eme
#
###################################################################################
use lib '/home/jmartijn/scripts/perllib';

#use lib '/local/one/LAUME/SOFTWARES/bin/perllib';
use lauralib ;
use Bio::SeqIO::fasta;
use Bio::SearchIO;
use Bio::DB::Taxonomy;
use Bio::DB::Taxonomy::flatfile;

use Scalar::Util qw(blessed);

$USAGE = "\n fasta2taxo_complete2.pl version 1.0 (Feb 8, 2016)\n\n [USAGE]   : fasta2taxo_complete2.pl [fasta file]\n\n" ;

# Print USAGE
unless( @ARGV )
{
    print $USAGE ;
    exit ;
}

$ali = $ARGV[0];
$out = $ali.".taxo";

### Open the fasta file
$searchio = Bio::SeqIO->new(-format => 'fasta', -file => $ali);

#print "Hey!\n";

$idx_dir = '/local/one/db/tax/';
($nodefile,$namesfile) = ('/local/one/db/tax/nodes.dmp','/local/one/db/tax/names.dmp'); #NCBI file containing taxonomy

#print "Hey hey!\n";

#my $obj = Bio::DB::Taxonomy::flatfile->new(-directory => $idx_dir,
#					  -nodesfile => $nodefile,
#					  -namesfile => $namesfile,
#					  -force     => 0);

$db = new Bio::DB::Taxonomy(-source => 'flatfile',
			    -nodesfile => $nodefile,
			    -namesfile => $namesfile,
			    -directory => $idx_dir); # Bioperl function parsing the NCBI taxonomy file


while( my $seqobj = $searchio->next_seq ) {
    
	# Reinitialize variables
	$id = "";
	$desc = "";
	$new_id = "";

	# retrieve variable for each sequence
	$id = $seqobj->display_id();
	$desc = $seqobj->description();
	$line = $id." ".$desc;

	$species = "";
	$genus = "";
	$complete_taxo = "";
	$switch_genus = 0;
	@taxo = ();

	if (($line =~ /YP/) || ($line =~ /NP/) || ($line =~ /XP/) || ($line =~ /gi/) || ($line =~ /\[/)) {
        
		@tab = split(/\]/, $line);
		$genus = $tab[-1];
		$annotation = $tab[0];        
		@tab = split(/\|/, $annotation);
		$annotation = $tab[-1];
		@tab = split(/\[/, $annotation);
		$annotation = $tab[0];
#		print $annotation,"\n";

		# retrieve genus and species name
		@tab = split(/\[/, $genus);
		$genus = $tab[1];
		($genus,$species) = split(/ /, $genus);
		$taxon = $genus." ".$species;

	#	print " --------> ", $taxon, "\n";

		# retrieve accession number
		@tab = split(/\|/, $line);
		$new_id = $tab[3];
		($new_id) = split(/\./,$new_id);

		# Retrieve classif
		@taxonids = $db->get_taxonids($taxon);
		$taxon_id = $taxonids[0];


		#       $leaf_classif = $db->get_taxon(-name => '$taxon');
		$leaf_classif = $db->get_taxon(-taxonid => $taxon_id);
		$obj_type = blessed($leaf_classif);

		if ( ($obj_type eq '') ){
			$leaf_classif = $db->get_taxon(-name => $genus);
			$obj_type = blessed($leaf_classif);
		}


		if ( ($line =~ /mixed culture/) || ($obj_type eq '') || ($genus =~ /synthetic/) || ($genus =~ / vector/) || ($genus =~ /uncultured/) || ($genus =~ /unidentified/) || ($genus eq 'unclassified') || ($genus eq 'Cloning') || ($genus eq 'Expression') || ($line =~ 'enrichment') || ($line =~ 'lasmid') || ($line =~ 'Shuttle') ) {
		$new_header = $genus."_".$species."@".$new_id;
		}
		else{
			
			$c =0;
				
			LOOP_A:until ($leaf_classif->rank eq 'superkingdom') {
				$leaf_classif = $db->ancestor($leaf_classif);
				if ($leaf_classif->scientific_name =~ /^\W/){
					#last LOOP_A;
					last;
				}
				
				$c++;
				push(@taxo, $leaf_classif->scientific_name);
			}
				
			$new_header = $genus."_".$species."@".$new_id;
		
			#$count_taxo = $c-1;
			$count_taxo = $c;
			$c_taxo = 1;
			
#			print "count taxo: ".$count_taxo."\n";
			
			if ($count_taxo >= 1 ){
				
				# We only want 4 taxonomical ranks
				if ($count_taxo > 10){
					$count_taxo = 10;
				}
				
#				print "c_taxo: ",$c_taxo,"\n";
				until ( $c_taxo == $count_taxo ){
					$new_header = $new_header."_".$taxo[-$c_taxo];
					$c_taxo++;
				}
			}
		}#ed else
	
	    $new_header = $new_header."_zzz_".$annotation;
	}
    
	else{
		$new_header = $line;
	}

	#### replace unwanted characters
	$id = $new_header;

	$id =~ s/:/_/g;
	$id =~ s/\"/_/g;
	$id =~ s/\'/_/g;
	$id =~ s/-/_/g;
	$id =~ s/;//g;
	$id =~ s/,//g;
	$id =~ s/\[//g;
	$id =~ s/\]//g;
	$id =~ s/\(//g;
	$id =~ s/\)//g;
	$id =~ s/\.//g;
	$id =~ s/ /_/g;
	$id =~ s/=/_/g;
	$id =~ s/\|/_/g;
	$id =~ s/\//_/g;
	$id =~ s/_+/_/g;
	$id =~ s/_$//g;

	print ">",$id,"\n";
	# sequence itself
	$sequence = $seqobj->seq();
	print $sequence,"\n";
	
}


close(OUT);

