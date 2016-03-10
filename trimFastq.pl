#!/usr/bin/perl

use strict;
use warnings;

my $usage = "\nUsage: $0 input.fq 5-prime-bp 3-prime-bp

     This script removes a specified number of bp from the 5' and 3' ends of each sequence in the fastq file\n\n";

# print "1: $ARGV[1]\n";
# print "2: $ARGV[2]\n";

defined $ARGV[2] or die $usage;

my $begTrim = $ARGV[1];
my $endTrim = $ARGV[2];

open( FAS, $ARGV[0] ) || die "Couldn't open $ARGV[0]: $!\n";

while( my $line = <FAS>){
  if( $line =~ /\A@/ ){
    print $line and next; # skip headers
  }elsif( $line =~ /\A\+/ ){
    print $line and next; # skip "+" line
  }else{
    my $len = length $line;
    print substr( $line, 0 + $begTrim, $len - $endTrim - $begTrim - 1 ), "\n";
  }
}

close FAS;

exit;
