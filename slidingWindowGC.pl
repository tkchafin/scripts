#!/usr/bin/perl 


use strict; 
use warnings; 
use Getopt::Long; 


our $input=""; 
our $win=100; 
our $inc=50;

parseArgs(); #Call subroutine to parse arguments... 

my $dna; 
my $name = "";  
my $header = ""; 
open ( FAS, "$input" ) || die "\nDerp: Can't open $input!\n$!\n"; 


#This block submits slidingWindowGC for each separate sequence in fasta file.
while (<FAS>) { 
    chomp $_; 
    if ( $_ =~ m/^>(\w+)/){ 
	$header = "$1";  #New sample name stored
	if ( $dna ) {
	    print "\n$name\n\n";
	    slidingWindowGC( $dna, 0 ) ;
	}     
	undef $dna;
	$name = $header; 
 
        }else{
	    $dna .= $_;  
	}
    }     

print "\n$name\n\n"; 
slidingWindowGC( $dna, 0 ); 



close FAS; 




#############################################SUBROUTINES###############################################

#Subroutine to parse command line arguments
sub parseArgs{

    my $usage = "\nUsage: $0 --input=file.fasta --window=100 --increment=50

    mandatory
       --file        -  FASTA file containing sequences; the first sequence in the file will be used
       --window      -  window length (default=100) 
       --increment   -  increment length; how far to shift each window (default=50) \n\n";


                 my $result = GetOptions
                         (
                                 'file=s'     => \$input,
				 'window=s'   => \$win,
                                 'increment=s'=> \$inc,
                               
                         );
             
	        $input eq "" and die $usage;  #Die if mandatory variables undefined
		$win==100 and print "\nWarning: No window length defined- using default of 100\n\n"; 
		$inc==50 and print "Warning: No increment length defined- using default of 50\n\n";
        
}
                                                             


#Recursive subroutine to perform sliding window through input DNA sequence

sub slidingWindowGC{ 
 
 
my $DATA = $_[0]; 
my $subseq;  
my $GC; 
my $start=$_[1];  # $start initialized at zero 

$subseq = substr ($DATA, $start, $win); #the "window"... 

$GC =()=$subseq =~ /G|C/gi; #Count up Gs and Cs 

print "$start\t$GC\n";  #print window start coordinate and GC content

$start+=$inc; #Increment start. The "sliding" part

#Check if $start is within length of the dna... Sets limit to recursive subroutine and keeps it from going crazy
    if ($start < length($DATA) ){ 
        slidingWindowGC( $DATA, $start);

    }
}






