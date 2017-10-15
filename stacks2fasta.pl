#!/usr/bin/perl 

use strict; 
use warnings;
use Getopt::Long;
use File::Path; 

our $input; 
our $workdir="";
our $catalog;
our $batch=1; 

parseArgs(); 
 
my $locus;
my $name; 
my @info; 
my %whitelist;
my $ID;
my $output="loci.$batch";
$workdir =~ /\S/ and $output = "$workdir/$output"; 

#Build list of loci containing SNPs (batch_#.catalog.snps.tsv output from STACKS) 
open ( CAT, $catalog ) || die "\nDerp: Can't open $catalog!\n\n";

while (<CAT>){ 
    @info = split /\t/, $_; 
    $ID = $info[2];

	#If locus ID is already in hash, then skip it 
    if (exists $whitelist{$ID}){ 
	 next;
    }else{
	$whitelist{$ID}=""; 
    }
}

close CAT;

#Parse STACKS output fasta file into loci, query each locus against whitelist
open ( IN, $input ) || die "\nDerp: Can't open $input!\n\n"; 

rmtree $output;
mkdir $output;
chdir $output;
 
while (<IN>){ 

$_ =~ m/CLocus_(\d+)_Sample_(\d+)/;
 
$locus = $1;
$name=$2;  

    if (exists $whitelist{$locus}){
	open ( OUT, ">>$locus.fasta");
	if ( $_ =~ />/ ){ 
	    print OUT ">$name\n";
	}else{
	    print OUT $_; 
	}
    }
}





##########################################SUBROUTINES#########################################

sub parseArgs{

	my $usage="\nstacks2fasta.pl takes the fasta output from STACKS and outputs a new fasta file for each locus containing variation, which are identified by querying the cstacks catalog

Usage: $0 --i /path/to/infile --w /path/to/workdir --c=/path/to/catalog

Mandatory Variables 
	-i, --input	-   path to input file (absolute path)
	-w, --workdir	-   path to working directory (new fasta files will be placed within /workdir/loci
	-c, --catalog	-   path to STACKS catalog 

Optional
	-b, --batch	-   Provide a batch number to append to output dir name [default=1]\n\n";

	my $result = GetOptions 
	( 
	'input|i=s'	=> \$input, 
	'workdir|w=s'	=> \$workdir, 
	'catalog|c=s'	=> \$catalog,
	'batch|b=i'	=> \$batch,
	); 

if ( $input eq "" ){ die "\nDerp: Input not specified!\n\n$usage"};

}

############################################################################################



