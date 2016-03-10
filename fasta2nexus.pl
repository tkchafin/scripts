#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Path;
use File::Basename;
use Parallel::ForkManager; 

# Declare variables

our @input;     
our $PROCS=1; 
#our $infiletype=1; 

parseArgs(); 
    
my ( $filepath, $dirpath ) = fileparse($input[0]);    
my $pm = Parallel::ForkManager->new($PROCS); 

#Iterate through files

@input = glob "@input";

foreach my $file ( @input ){ 

  $pm->start and next; #Begin one daughter process, start the next.
	
#Initialize variables within each daughter process
    my @data;
    my @names;
    my $taxa = 0;
    my $name; 
    my @fasta;
    my @loci; 
    my $nchar; 
    my $line; 
 
    open ( FILE, "$file" ) || die "Do I need to call the short bus?\nCan't open $file: $!\n";
	while ( <FILE> ){
			chomp;
				 
			if( $_ =~ /^\>/ ){
				$taxa++; 
				$_ =~ /^\>(\S+)/;
				$name = "$1";
				push @names, "$name";
			}elsif( $_ =~ /^\s*$/ ){
				next;
			}elsif( $_ =~ /^\s*#/ ){
				next;
			}else{
				push @data, $_ ;
				$nchar = length;
			}
    	}
    close FILE;

    #Capture taxa name to use as identifier
	my $filepath = fileparse("$file");
	$filepath =~ /(\w+)\.\w/;
	my $ID = $1;

    open( OUT, '>', "$dirpath$ID.nex" ) || die "Do I need to call the short bus?\nCan't write to $ID.nex\n";		
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
  $pm->finish;   #Indicate that daughter process is complete
}


$pm->wait_all_children;   #Wait for all daughter processes before exiting program.

exit;
###########################SUBROUTINES###################################

sub parseArgs{ 
	#Message to print if mandatory variables not declared
	my $usage ="\nUsage: $0 --i /path/to/input/directory/*.fasta
Mandatory 
	-i, --input	-  path to the input files in fasta format 
	-p, --procs	-  Number of processes for multithreading [default=1]
\n";

	my $options = GetOptions 
		( 
		'input|i=s{1,}'		=>	\@input,
		'procs|p=i'		=>	\$PROCS,
		);
		
	@input or die "\n\nDo I need to call the short bus?: Input not specified!\n\n$usage\n"; 
}

#########################################################################
