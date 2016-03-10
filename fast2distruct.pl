#!/usr/bin/perl 

use strict; 
use warnings;
use Getopt::Long;

my $struct;  
my $popmap ;
my $meanQ  ;
my $out = "./distruct";
my $inline = 0;
my $substrGuess = 0; 
my $help = 0; 

parseArgs(); 

my $popq = $out . ".popq"; 
my $indivq = $out . ".indivq"; 

#Capture individual order from structure file
open (my $str, $struct) || die "\nDerp: Couldn't open $struct!\n\n";

my @indOrder;
my @row; 
my $index = 0; 
my $skipLine = 0; 
while (<$str>){
  #If 2 line structure format, skip every other line
  if ($inline == 0){
    if ($skipLine == 1){
      $skipLine = 0;
      next; 
    }
  }
  chomp; 
  @row = split /\t/, $_;
  s/\s+//g; 
  next unless length; 
  $row[0] =~ s/\s+//g;
  $indOrder[$index] = $row[0]; 
  $index++;
  $skipLine = 1;  
}
close $str; 

#Capture population identifiers 
open (my $pops, $popmap) || die "\nDerp: Couldn't open $popmap!\n\n"; 

my %popHash;
my %substrHash; 
my %seen; #Hash lookup table 
my $popcount = 1; 
my $sub;
my $subCount = 1;  
while (<$pops>){
  chomp; 
  @row = split /\t/, $_; 
  $row[0] =~ s/\s+//g;
  $row[1] =~ s/\s+//g; 
  s/\s+//g; 
  next unless length; 
  $row[0] = uc($row[0]); 
  $row[0] =~ /(\d+[A-Za-z]+)/; 
  $sub = $1;  
  if ($seen{$row[1]}){ 
    $popHash{$row[0]} = $seen{$row[1]};
    if ($substrGuess == 0){
      if ($substrHash{$sub}){ 
        $substrHash{$sub} != $seen{$row[1]} and print "Warning:Found more than one population identifier for the same population substring " . $sub ."!\n"; 
      }else{
        #print "2Setting " . $sub ." to " . $seen{$row[1]} . "\n";
        $substrHash{$sub} = $seen{$row[1]}; 
      }
    }
  }else{ 
    $seen{$row[1]} = $popcount;
    $popHash{$row[0]} = $popcount; 
    $substrGuess == 0 and $substrHash{$sub} = $popcount; 
    #$substrGuess == 0 and print "1Setting " . $sub ." to " . $popcount . "\n";
    $popcount++; 
  }
    
} 
close $pops; 

#Parse meanQ file and write indivq file and make calculations for popq
open (my $results, $meanQ) || die "\nDerp: Couldn't open $meanQ!\n\n"; 
open (my $iq, ">$indivq") || die "\nDerp: Couldn't open $indivq!\n\n";

my @asn;
$index = 0;  
my $popID; 
my %popq; 

while (<$results>){ 
  chomp;
  @row = split /\s+/, $_;
  s/\s+//g; 
  next unless length; 
  $asn[$index] = [@row];  

  if ($popHash{uc($indOrder[$index])}){ 
    $popID = $popHash{uc($indOrder[$index])}; 
  }elsif ($substrGuess == 0){
    print "\nWarning: Individual ". $indOrder[$index]; 
    print " not found in popmap! Trying to guess correct population identifier... \n"; 
    $indOrder[$index] =~ /(\d+[A-Za-z]+)/; 
    $sub = uc($1); 
    if ($substrHash{$sub}){ 
      $popID = $substrHash{$sub};
      print "Assigning " . $sub ." to population " . $substrHash{$sub} . "\n";
    }else{
      print "Population substring " . $sub . " not found. Setting popID to ".$popcount .".\n";
      $popID = $popcount;
      $popHash{$indOrder[$index]} = $popcount; 
      $substrHash{$sub} = $popcount; 
      $popcount++;  
    }
  }elsif ($substrGuess == 1){ 
    
    print "\nWarning: Individual " . $indOrder[$index]; 
    print " not found in popmap! Setting popID to ". $popcount . ".\n";
    $popID = $popcount;
    $popHash{$indOrder[$index]} = $popcount; 
    $popcount++;  
  }
  print $iq "  " . $indOrder[$index] .  "      " . $index . "   (0)     "; 
  print $iq $popID. "  : "; 
  printf $iq "%.4f ", $_ for @row; 
  print $iq "\n"; 
  
  #If popID already in popq table, then add to it
  if ($popq{$popID}){ 
    for (my $i =0; $i <= $#row; $i++){ 
      $popq{$popID}[0][$i] += $row[$i];  
    }
    $popq{$popID}[1]++; 
  }else{ 
    $popq{$popID}[0] = [@row]; 
    $popq{$popID}[1] = 1; 
  }

  $index++; 
}
close $meanQ; 
close $iq; 

#Process popq table and print popq file 
open (my $pq, ">$popq") || die "Derp: Couldn't open $popq!\n";
my $total; 
foreach my $key ( sort {$a<=>$b} keys %popq){ 
   print $pq $key . ":   ";
   $total = $popq{$key}[1]; 
   for (my $i=0; $i<@{$popq{$key}[0]}; $i++){ 
      printf $pq("%.4f",($popq{$key}[0][$i]/$total));
      print $pq " ";
   }
   print $pq $total . "\n"; 
}

close $pq; 
exit;

#########################################################################################

sub parseArgs{

my $message = 
"This script converts from the output of fastStructure to the input required for standard Distruct. It requires the structure file output by pyRAD (which was used for the analyses) and a population map in the style SMM required for Astral pipeline. It will use that pop map to determine a priori groupings, for building the popq files. I might add the ability to just pull these from the structure file later, but the pyRAD str doesn't have this so that's why I didn't do that yet. 

If you have problems running the script let me know. It hasn't really been tested fully, and I threw it together quickly. 

Options 

	-i	- Input fastStructure meanQ file 
	-s	- Structure file from pyRAD which was given to fastStructure 
	-p	- Population map, tab-delimted (e.g. 8HBC001 \t cypha or 8HBC001 \t 1) 
	-o	- Output prefix (e.g. ./k3) 
	-e	- Bool, toggle to turn off population estimation based on prefix 
		  e.g. program will guess 9WRW002 goes in same population as 9WRW001
		  if 9WRW002 is missing from pop map. I did this because my pop map was 
 		  missing samples, and I didn't want to go back and fine all of them to add in 

"; 

	my $result = GetOptions
	( 
	'i=s'	=> \$meanQ, 
	's=s'	=> \$struct, 
	'p=s'	=> \$popmap, 
	'o=s'	=> \$out, 
	'h!'	=> \$help,
	'e!'	=> \$substrGuess 
	);
$meanQ or die $message; 
$struct or die $message; 
$popmap or die $message; 
$help == 1 and die $message; 

}


 
