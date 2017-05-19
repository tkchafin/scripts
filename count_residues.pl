#!/usr/bin/perl 

#Author: Tyler K. Chafin
#Script counts a given residue in provided amino acid alignment

use strict;
use warnings;

my $file = $ARGV[0];
my $res = $ARGV[1];

open (FASTA, "$file") || die "Cannot open file $1: $!\n";
my %outhash;
my $temp;
while (<FASTA>){
	chomp $_;
	if ($_ =~ "^\>"){
		$temp = $_;
		$outhash{$temp} = "";
		next;
	}else{
		$outhash{$temp} .= $_;
	} 

}
close FASTA;

open (OUT, ">out.tsv") || die "Cannot open out.tsv: $!\n";
for my $key (keys %outhash){
	$key =~ m/\>.*\|.*\|(.*?)\s+.*/;
	my $match = $1;
	my $count = () = $outhash{$key} =~ /$res/gi;
	print OUT $match, "\t", $count ,"\n";
}
close OUT;

exit;
