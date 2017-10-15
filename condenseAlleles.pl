#!/usr/bin/perl 

use strict; 
use warnings;
use Getopt::Long;
use File::Path; 
use File::Basename; 

our @input; 

#Call subroutine to parse command-line arguments
parseArgs(); 
 


#iterate through input files 
@input = glob "@input";
foreach my $file ( @input ){ 

my $name = "";
my %data;

#Open each file
    open ( FAS, "$file" ) || die "Derp: Cannot open $file!";
 
    while (<FAS>){ 

	chomp $_; 
	
	if ($_ =~ /\>/ ){ 
	    $_ =~ />(\S+)/;
 	    $name = "$1";
	}else{ 
	    push @{$data{$name}}, $_;
	}
    }
	consense( \%data, $file ); 
    
}

close FAS; 

###########################################SUBROUTINES#########################################

sub parseArgs{

	my $usage="\ncondenseAlleles.pl takes FASTA files with individual alleles encoded separately, and generates a consensus sequence for each individual

Usage: $0 --i /path/to/*.fasta

Mandatory Variables 
	-i, --input	-   path to input files in FASTA format\n\n"; 

	my $result = GetOptions 
	( 
	'input|s=s{1,}'	=> \@input,
	); 

@input or die "\nDerp: Input not specified!\n\n$usage";

}

############################################################################################

#This subroutine takes a hash of arrays, where the array stores the alleles for each sample (key in the hash), and spits out a consensus sequence for any heterozygote. 
#
#WARNING: This script may contain excessive warnings!
#
#WARNING: OVERWRITES THE ORIGINAL FILE!!
#
#ANOTHER WARNING: Not written to accomodate polyploids, or paralogous loci. I will be implementing this as a part of post-processing the output of Stacks, which should have removed any loci with individuals having more than 2 alleles (as these indicate presence of paralogs, which I don't want). 
#
#YET ANOTHER ANOTHER WARNING: This script also assumes that alleles are of the same length, and pre-aligned. 

sub consense{

my $datref = $_[0]; 
my $file = $_[1];




#Empty current file to start rewriting (all info stored in %info hash now) 
open ( OUT, ">$file" ) || die "Derp: Oh noes! I'm in the subroutine and cannot open $file!"; 

foreach my $key ( sort {$a <=> $b } ( keys %$datref ) ) { 
    
    print OUT "\>$key\n";  #Print FASTA header 

    if ( exists $$datref{$key}->[1] ) { #Check if sequence has multiple alleles
	#Load sequences into arrays for comparison
	my @allele1 = split //, $$datref{$key}->[0];
	my @allele2 = split //, $$datref{$key}->[1]; 

	for ( my $i=0; $i <= length $$datref{$key}->[0]; $i++ ){
		#This is messy and sucks. Fix later with Bio::AlignIO and consensus_iupac function
	    if ( uc ( $allele1[$i] ) eq uc ( $allele2[$i] ) ) { 
		uc ( $allele1[$i] ) eq "A" and print OUT "A"; 
		uc ( $allele1[$i] ) eq "G" and print OUT "G"; 
		uc ( $allele1[$i] ) eq "T" and print OUT "T"; 
		uc ( $allele1[$i] ) eq "C" and print OUT "C"; 
	    }else{ 
		if ( uc $allele1[$i] eq "A" ){	
		    uc $allele2[$i] eq "G" and print OUT "R"; 
		    uc $allele2[$i] eq "C" and print OUT "M"; 
		    uc $allele2[$i] eq "T" and print OUT "W"; 
		}elsif ( uc $allele1[$i] eq "G" ){ 
		    uc $allele2[$i] eq "A" and print OUT "R"; 
		    uc $allele2[$i] eq "C" and print OUT "S"; 
		    uc $allele2[$i] eq "T" and print OUT "K"; 
		}elsif ( uc $allele1[$i] eq "T" ){ 
		    uc $allele2[$i] eq "A" and print OUT "W";  
		    uc $allele2[$i] eq "C" and print OUT "Y";
		    uc $allele2[$i] eq "G" and print OUT "K";
		}elsif ( uc $allele1[$i] eq "C" ){
		    uc $allele2[$i] eq "A" and print OUT "M";  
		    uc $allele2[$i] eq "G" and print OUT "S"; 
		    uc $allele2[$i] eq "T" and print OUT "Y";
		}
	    }
	} 
	print OUT "\n"; 
    }else{ 
	print OUT $$datref{$key}->[0],"\n";
    }    
}
#close OUT;
}

