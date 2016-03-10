#!/usr/bin/perl 

#Tyler K. Chafin; 14-Dec-15 
#tkchafin@uark.edu 

use strict; 
use warnings; 

my $usage = "
This script functions to collapse aligned sequences in FASTA format into haplotypes, and sort halpotypes by frequency. 

Author: Tyler K. Chafin - tkchafin\@uark.edu
Last Modified: 14-Dec-15

Usage: $0 </path/to/.fasta> 

"; 
my $file; 
if (defined $ARGV[0]){
  $file = $ARGV[0]; 
  print "Input: $file\n"; 
}else{ 
  die $usage; 
}

my %contents; 
my %freq; 

open (INPUT, $file) || die "Cannot open $file: $!\n\n"; 
while (<INPUT>){ 
  chomp;   
  if ($_ =~ /^\s*$/){
    next;  
  }elsif ($_ =~ /^>/){ 
    next; 
  }else{
    if ($contents{uc $_}){
      $contents{uc $_}++; 
    }else{ 
      $contents{uc $_} = 1; 
    }
  }   
}
close INPUT;
open (OUT, ">sorted.fasta") || die "Could not open output file: $!\n"; 
print "Output: sorted.fasta\n";
my $count=1; 
foreach my $key (sort {$contents{$b} <=> $contents{$a}}keys %contents){ 
  print OUT ">H$count\n";
  print OUT "$key\n";
  $count++; 
}
close OUT;

exit; 
