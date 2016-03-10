#!/usr/bin/perl

# gives the non-gap-character-length of each sequence

use warnings;
use strict;

die "usage: $0 fastafile\n" unless $#ARGV == 0;

open A, shift;

my ($id, @ids, %seq, $total);
while (<A>) {
  chomp;
  if (/^>(.*)/) {
    $id = $1;
    push @ids, $id;
  } else {
    $seq{$id} .= $_;
  }
}

my (%group, $group, $ar, $seq, $len);
foreach $id (@ids) {
  $seq = $seq{$id};
  $seq =~ s/[\s*-]+//g;
  $len = length $seq;
  $total += $len;
  print "$id $len\n";
#  print "\t$len\n";
#  print ">$id ; length $len\n$seq{$id}\n";
}

print "$total total sequence length\n";
