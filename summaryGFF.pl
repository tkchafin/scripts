#!/usr/bin/perl 


use strict; 
use warnings; 
use Getopt::Long; 


our $gff=""; 
our $genome=""; 

parseArgs(); #Call subroutine to parse arguments... 

my @line; 
my $dna; 
my $subseq; 
my $total;
my $GC;  
my $element;
my $add; 
my @info;  
my %summary; 

	



#Call subroutine for each type of element...



summaryGFF(); 







################################SUBROUTINES######################################

#Subroutine to parse command line arguments
sub parseArgs{

    my $usage = "\nUsage: $0 --genome=whole_genome.fasta --gff=annotations.gff

    mandatory
       --genome      -  FASTA file containing sequences to parse
       --gff         -  GFF file containing gene annotations \n\n";


                 my $result = GetOptions
                         (
                                 'genome=s'  => \$genome,
                                 'gff=s'     => \$gff,
                               
                         );
             
	        $genome ne "" || die $usage;  #Die if mandatory variables undefined
		$gff ne "" || die $usage; 
        
}
                                                             

#Subroutine to parse gff and genome for particular type of element

sub summaryGFF{


undef @line; 
undef $dna; 
 
    
    open ( GENOME, "$genome") || die "Derp: Can't open file $genome!";

	while (<GENOME>){ 
	    $_ ne /^>/ and $dna .= $_; 
	};

    close GENOME;
 

    open ( GFF, "$gff" ) || die "Derp: Can't open file $gff!"; 

	foreach ( <GFF> ){ 
		@line = split /\t/, $_;
			#print "$line[2]\n"; 
	            $GC=0; 
		    $subseq = substr ( $dna, $line[3]-1, $line[5] ); 
                    $add =()=$subseq =~ /G/gi; 
		    $GC += $add; 
		    $add =()=$subseq =~ /C/gi;
		    $GC += $add; 
		        #print "$GC\n";

#If element is already in hash, then alter values in the arrays by following ref in hash value...
		    if ( exists $summary{$line[2]} ){ 
		            #print "$line[2]\n";
			$summary{$line[2]}->[0] += $line[5];
			   #print "$summary{$line[2]}\n"; 
			$summary{$line[2]}->[1] += $GC;
		    }else{ 
		
 #Create array containing length and GC content, then assign array ref to hash key for that element	
			my @info=($line[5], $GC); 
			$summary{$line[2]} = \@info; 
		    } 	 

	}
	
	foreach my $key ( keys %summary ){ 
	    print "$key \t$summary{$key}->[0]   ";
	    printf( "(%.1f%%) \t", $summary{$key}->[0] / length($dna) * 100); 
	    printf( "%.2f \n", $summary{$key}->[1] / $summary{$key}->[0] * 100); 
        }
}	

close GFF; 
exit; 



