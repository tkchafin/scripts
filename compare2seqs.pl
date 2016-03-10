#!/usr/bin/perl

use strict; 
use warnings; 

my $tax1="Carica"; 
my $tax2="Vitis"; 
my $input="cox1.fasta";

#Open input fasta
open (FAS, "$input" ) || die "\n\nI pity the fool that can't open their fasta files\n\n$!\n\n";

my $count=0; 
my @dna1; 
my @dna2;


#Set first and second sequences to exploded arrays
while ( <FAS> ){ 
    if ($_ !~ /^>/)  { 
	$count++;  
	$count==1 and @dna1 = split //, "$_";
	$count==2 and @dna2 = split //, "$_"; 
    }
}
		  
#print "\n\n@dna1\n\n@dna2\n\n"; 


if (length(@dna1) ne length(@dna2)){
    print "\nWarning:Sequences to compare are of different length. Check your alignment.\n\n";
}

#Print header for table
print "\nPosition $tax1 $tax2\n"; 

my $samelen=0;
my $difflen=0;
my $total=0; 
for ( my $i=0; $i <= $#dna1; $i++){ 
    if ($dna1[$i] eq $dna2[$i]){ 
	$samelen++;
	$total++; 
    }else{ 
	$difflen++;
	$total++;  
	print $i+1, "\t$dna1[$i]\t$dna2[$i]\n"; 
    } 
} 

print "\nnumber of identical sites: $samelen\n"; 
print "number of different sites: $difflen\n"; 
print "percent difference: "; 
printf( "%.2f \n\n", $difflen / $total * 100 ); 




