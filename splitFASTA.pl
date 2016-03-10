#!/usr/bin/perl

#
#Script by Tyler K. Chafin
#Last Modified: 6 May 2015
#Added: Capability to split file into user-defined number of parts
#


use strict; 
use warnings;
use Getopt::Long;

our $pattern=">";
our $infile="";
our $suffix="fasta";
our $breaks;  
parseArgs();



open ( INFILE, "$infile" ) ; 

my $n=0; 
my $matches=0; 
if ($breaks){
#print "Breaks = $breaks\n"; 


    #Count num of pattern matches
    while (<INFILE>){ 
	$_ =~ "$pattern" and $matches++; 
    }  	
    my $num_lines = int($matches/$breaks); 
    my $count=0;
 
#print "Num_lines = $num_lines\n"; 
#print "matches = $matches\n"; 
    #Foreach part, read data and write to appropriate outfile 
     
    seek(INFILE,0,0); #Reset reading position in fh
    
    while (<INFILE>){ 
 
	if ($_ =~ "$pattern"){ 
	    $n++;
	    if ($count == 0){ 
		$count++;
		open (OUTFILE, "> $count.$suffix") || die $!;  
	     
	    }
	    if ($count >= $breaks){ 
		print OUTFILE "$_"; 
	 
	    }else{ 
	    	if($n<=($num_lines*$count)){ 	
		    print OUTFILE "$_"; 
	    	}else{ 
		    close OUTFILE; 
		    $count++; 
		    open (OUTFILE, "> $count.$suffix") || die $!;
		
		    print OUTFILE "$_"; 
	    	}
	    } 
	}else{ 
	    print OUTFILE "$_"; 
	} 
   }		

#If num breaks not defined: break for each contig
}else{
    while (<INFILE>) { 

    	if ( $_ =~ "$pattern"  ){
            $n++; 
	    open ( OUTFILE, "> $n.$suffix" ) || die $!;  	
            print OUTFILE "$_"; 	
    	}else{
	    print OUTFILE "$_";
        } 
    }
}


close INFILE;
close OUTFILE; 

exit;
###################################################

sub parseArgs{ 

    my $usage = "\nUsage: $0 --file=whole_genome.fasta --pattern=\> --suffix=fasta

Author: Tyler K. Chafin - tkchafin\@uark.edu
Last Modified: 6 May 2015

Purpose of script is to break a given FASTA file into a user-defined number of portions, or into separate files per FASTA header and associated sequence. 


mandatory
   --file      -  File to break up

optional 
   --breaks    -  Break file into n pieces [default is one file for each contig]
   --pattern   -  Pattern to use to divide file [default=>] 
   ---suffix   -  Suffix to use when naming daughter files [default=fasta]\n\n";


    my $result = GetOptions 
	(
	'f|file=s'	=> \$infile, 
	'p|pattern=s'	=> \$pattern, 
	's|suffix=s'	=> \$suffix, 
	'b|breaks=i'	=> \$breaks, 
	); 
    
    if ( $infile eq "" ){ die $usage}; 
}

