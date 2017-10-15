#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Path;

# Declare variables
    
    my $line; 
    my @loci;
    my $workdir=""; 
    my @fasta;
    my $i = 1;
    my $input;
    my $batch;
	 
parseArgs();

my $output="loci";
$batch and $output.=$batch;

# open file and read it in
open( LOCI, $input ) || die "Can't open $input: $!\n";
# Make the loci directory to put fasta files from pyrad2fasta subroutine
chdir $workdir;
rmtree $output;
    mkdir $output;
    chdir $output;
    
while ( $line = <LOCI> ){

    
    if( $line =~ /^\/\// ){
        if( $line =~ /\*|\-/ ){
            pyrad2fasta( @loci, $i );
			undef( @loci );
			$i++;

        }else{
            undef( @loci );
        }	
    }else{
        push @loci, $line;
		
    }
    

}
close LOCI;
exit;

###########################SUBROUTINES###################################

sub parseArgs{ 
	#Message to print if mandatory variables not declared
	my $usage ="\npyrad2fasta.pl takes the custom .loci output from pyRAD and creates a new FASTA file for each locus containing at least 1 SNP.

Usage: $0 --i /path/to/*.loci --w /path/to/workdir   

Mandatory 
	-i, --input	-  path to the input file (*.loci from pyRAD)  
	-w, --workdir	-  path to working directoy (new fasta files will be placed within /workdir/loci 

Optional
	-b, --batch	- Provide a batch number

\n";

	my $options = GetOptions 
		( 
		'input|i=s'		=>	\$input,
		'workdir|w=s'		=> 	\$workdir,
		'batch|b=i'		=>	\$batch,
		);
		
	$input or die "\n\nDo I need to call the short bus?: Input not specified!\n\n$usage\n"; 
	if ( $workdir eq ""){die "\nDerp: Working directory not specified!\n\n"}; 
}

#########################################################################

sub pyrad2fasta{

		
	# split at whitespace
    for my $element ( @loci ){
		open( OUT, '>>', "$i.fasta" ) || die "Do I need to call the short bus?  Can't write to $i.fasta: $!\n\n";
		 @fasta = split( /\s+/, $element );
		 print OUT $fasta[0], "\n";
		 print OUT $fasta[1], "\n";
	
	}
#     Print the loci in FASTA format		
		
		   
}
