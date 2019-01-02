#!/usr/bin/perl

#Modified from  fasta2nexus.pl written by BTM
#-TKC

use strict;
use warnings;
use Getopt::Long;
use File::Path;
use File::Basename;
# Declare variables

our $input;      
#our $infiletype=1; 

parseArgs(); 
         

#Initialize variables within each daughter process
    my @data;
    my @names;
    my $taxa = 0;
    my $name; 
    my @fasta;
    my @loci; 
    my $nchar; 
    my $line=0; 
    my @linedata; 
 
    open ( FILE, "$input" ) || die "Error\nCan't open $input: $!\n";
	while ( <FILE> ){
	    chomp;
	    $line++; 
	    @linedata = split /\s+/, $_; 
	    s/\s+//g; 
	    length($_) or next;
	    $line == 1 and next;  
	    $taxa++; 
	    $name = $linedata[0];
	    push @names, "$name"; 
	    push @data, $linedata[1];
	    if ($nchar){ 
		length($linedata[1]) != $nchar and print "Error: Line beginning with $name has a different sequence length.\n"; 
	    }else{ 
		$nchar = length($linedata[1]);
	    }
    	}
    close FILE;

    #Capture to use as identifier
	my ($filepath, $dirpath) = fileparse("$input");
	$filepath =~ /(\w+)\.\w/;
	my $ID = $1;

    open( OUT, '>', "$dirpath$ID.nex" ) || die "Error\nCan't write to $ID.nex\n";		
         print OUT "#NEXUS\n\n";    
         print OUT "BEGIN DATA;
DIMENSIONS NTAX=$taxa NCHAR=$nchar;
FORMAT DATATYPE=DNA MISSING=? GAP=- ;

MATRIX\n";

	for ( my $i = 0; $i<scalar @names; $i++ ){	
		print OUT "$names[$i]\t$data[$i]\n";
	}
    print OUT ";\n";

    print OUT "END;\n\n";
    
    close OUT; 	


exit;
###########################SUBROUTINES###################################

sub parseArgs{ 
	#Message to print if mandatory variables not declared
	my $usage ="\nUsage: $0 --i /path/to/input/directory/*.phylip
Mandatory 
	-i, --input	-  path to the input files in phylip format
\n";

	my $options = GetOptions 
		( 
		'input|i=s{1,}'		=>	\$input
		);
		
	$input or die "\n\nError: Input not specified!\n\n$usage\n"; 
}

#########################################################################
