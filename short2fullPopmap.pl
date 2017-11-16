#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

my %opts;
getopts('i:s:c:h', \%opts);

if ($opts{h}){
	&help;
	die "Exiting because help menu was called.\n\n"
}

my ($names, $short, $num, $out) = &parseArgs(\%opts);

open (my $fhs, $short) or die "Can't open short\n";

my %hash;
while (my $row = <$fhs>){
	chomp $row;
	my @arr = split "\t", $row;
	if  (!exists $hash{$arr[0]}){
		$hash{$arr[0]} = $arr[1];
		#print $arr[0]," is ", $hash{$arr[0]}, "\n";
	}
}
close $fhs;


open (my $fhn, $names) or die "Can't open names\n";
open (my $outfh, ">$out") or die "Can't open output file for writing\n";

while (my $name = <$fhn>){
	chomp $name;
	my $n = substr $name, 0, $num;
	if (exists $hash{$n}){
		print $outfh $name, "\t", $hash{$n}, "\n";
	}else{
		print "$name ($n) doesn't match anything", "\n";
	}
}
close $fhn;
close $out;

exit;

###############################################################################
################################ Subroutines ##################################
###############################################################################

# subroutine to print help
sub help{

  print "\nLazy script to create full popmap from a prefix popmap\n\n";
  print "Program Options:\n";

  print "\t-i:\tText file with list of sample names\n";
  print "\t-s:\tTab-delimited prefix names\n";
	print "\t-s:\tNumber of characters used for prefix\n";
	print "\t-o:\tOutput file name\n";
	print "\t-h:\tBoolean. Calls help menu.\n\n";

}


# subroutine to parse the command line options
sub parseArgs{

  my( $params ) =  @_;
  my %opts = %$params;

	my $names = $opts{i} or die "File with sample names not given\n";
  my $short = $opts{s} or die "File with prefix popmap not given\n";
	my $num = $opts{c} or die "Number of characters not given\n";
	my $out = $opts{o} || "output.popmap";


  return( $names, $short, $num, $out);

}
