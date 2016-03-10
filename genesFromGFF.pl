#!/usr/bin/perl 


use strict; 
use warnings; 
use Getopt::Long; 


our $gff=""; 
our $genome=""; 

parseArgs(); #Call subroutine to parse arguments... 

my @line; 
my $dna; 
my $element;
my @info;   
	




#Call subroutine- recognizes each element type in gff and provides total length and GC content for each

seqsFromGFF( "CDS" ); 
seqsFromGFF( "rRNA" );  
seqsFromGFF( "tRNA" ); 


exit; 

###############################SUBROUTINES######################################

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

sub seqsFromGFF{

my $type = $_[0]; 
my %genes; 
my $subseq; 
my $name; 
my $exon; 


undef @line; 
undef $dna; 
 
    
    open ( GENOME, "$genome") || die "Derp: Can't open file $genome!";

	while (<GENOME>){ 
	    if ($_ !~ />/){ 
		chomp $_; 
		$dna .= $_; 
	    }
	}

    close GENOME;
 

    open ( GFF, "$gff" ) || die "Derp: Can't open file $gff!"; 

	foreach ( <GFF> ){ 
            @line = split /\t/, $_;
			#print "$line[2]\n"; 
	        if ( uc $line[2]  eq uc $type ){ 
		    $subseq = substr ( $dna, $line[3]-1, $line[5] ); 
                    
		    @info = split /\s/, $line[8]; 
		    $name = $info[1]; 
		    
                    #Reverse complement if on opposite strand
		    $line[6] =~ "-" and $subseq = revcom( $subseq ); 


		    if ( uc $info[2] eq uc "exon"){ 
			$exon = $info[3]; 
		    }else{
			$exon = 1; 
		    }

		    		 


#If element is already in hash, then alter values in the arrays by following ref in hash value...
		    if ( exists $genes{$name} ){  
			$genes{$name}->[$exon-1] = $subseq;
		    }else{ 
		
 #Create array containing length and GC content, then assign array ref to hash key for that element	
			my @seqs=(); 
			$genes{$name} = \@seqs;
		   	$genes{$name}->[0] = "";
			$genes{$name}->[$exon-1] = $subseq;  
		    } 	 
		}
	}

print "\nSequences for element type \"$type\": \n\n";	

foreach my $key ( keys %genes ){ 
    print "\>$key\n";
    print "@{$genes{$key}}"."\n";
    
}
        
close GFF; 
}	
 
###################################################################################################

sub revcom { 

my $DNA = reverse ( $_[0] ); 

$DNA =~ tr/ACGTacgt/TGCAtgca/; 

return $DNA; 

}
