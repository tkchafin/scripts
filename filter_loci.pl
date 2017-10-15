#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long; 
use File::Path; 
use File::Copy; 
use File::Basename;

my @DIR; 
my $cutoff;
my $N; 
my $blacklist=0;
my @contents;

parseArgs();  



#Use File::Basename to capture some info
my ($filepath, $dirpath) = fileparse( $DIR[0] );

#If user toggled "blacklist loci" on, then some configuration...
$blacklist =~ "1" and my $bname = "blacklist";
$blacklist =~ "1" and rmtree "$dirpath/$bname";
$blacklist =~ "1" and mkdir "$dirpath/$bname";
 

#Iterate through files

@DIR = glob "@DIR"; 
foreach my $file (@DIR){ 

$N=0; 
@contents="";

open (FILE, "$file"); 
    while (<FILE>){ 
	chomp $_;
	push @contents, $_; 
    }
	    $N +=()= "@contents" =~ /\>/g;
	    	    
	    if ( $blacklist eq "1"){ 
	        if ( $N < $cutoff ){ 
		    move("$file","$dirpath/$bname/");  
	        }
	    }else{
	
	        if ( $N < $cutoff ){ 
		    unlink "$file"; 
	        } 
            }
close FILE; 
}


##############################################SUBROUTINES###########################################

sub parseArgs{ 

	my $usage="\nfilter_loci.pl takes a directory full of fasta files, each representing a single locus, and deletes any loci with insufficient coverage across samples, using a user-specified cut off value.                  

Usage: $0 --i=/path/to/*.fasta --x=# [--b] 

Mandatory Variables 
	-i, --input		-   Path to fasta files 
	-x, --cutoff		-   Integer indicating the minimum number of samples to retain a locus
Options
	-b, --blacklist		-   Retain dropped loci in a blacklisted_loci directory
";

	my $results = GetOptions 
	( 
	'input|i=s{1,}'	=> \@DIR,
	'cutoff|x=i'	=> \$cutoff,
	'blacklist|b!'	=> \$blacklist,
	); 

@DIR or die "\nDerp: Input directory not defined!\n\n$usage";
$cutoff or die "\nDerp: Minimum coverage required to retain a locus not defined!\n\n$usage";

}
