#! /usr/bin/perl

# By Tyler K. Chafin
# Contact: tkchafin@uark.edu

use strict; 
use warnings; 
use Getopt::Std; 

#Die if no arguments given
if( scalar( @ARGV ) == 0 ){
  &help;
  die "No options given.\n\n";
}

#Parse arguments
my %opts;
getopts( 'i:o:hm:n:x:', \%opts );

# kill if help option is true
if( $opts{h} ){
  &help;
  die "Printing help menu.\n\n";
}

#get options 
my ($in, $out, $min, $max, $cap) = &parseArgs(\%opts); 

open (FASTA, "$in") || die "Could not open file $in: $!\n";
print "\nReading input file <$in>...\n";
open (OUT, ">$out") || die "Could not open file for output ($out) : $!\n";
print "Writing output to <$out>...\n";
my $base;
my $count;
my $num = 0;
while (<FASTA>){
	chomp;
	if ($_ =~ /^\>/){ #If header line
		my @line = split(/-/, $_);
		$line[0] =~ s/\>//g;
		$base = $line[0];
		$count = $line[1];
		$num == 1 and die "Error: Header line \"$_\" immediately follows another header line.\n";
		$num = 1;
		next;
	}else{
		$num == 2 and die "Error: Sequence line \"$_\" immediately follows another sequence line.\n";
		$num = 2; 
		
		if ($count < $min && $min != 0){
			print "Skipping <$base>: Depth (<$count>) is below minimum <$min>!\n";
			undef($count);
			undef($base);
			next;
		}elsif ($count > $cap && $cap != 0){
			print "Skipping <$base>: Depth (<$count>) is above maximum <$cap>!\n";
			undef($count);
			undef($base);
			next;
		}elsif ($max != 0){
			$count > $max and $count = $max;
		}
		for (my $i=1; $i <= $count; $i++){
			print OUT ">" . $base . "-" . $i . "\n";
			print OUT $_ . "\n";
		}
		undef($count);
		undef($base);
	}
}
print "Done!\n\n";
close FASTA;
close OUT;
exit;

 ########################### SUBROUTINES ###############################

 sub help{
	 
	print "\nThis perl script is written by Tyler K. Chafin - tkchafin\@uark.edu\n";
	print "\nInput should be a FASTA file of collapsed read clusters where -# at the end of the FASTA header for each sequence indicates the stack depth for the cluster.\n";
	print "\nNOTE: Stack depth counts start at 1.\n";
	print "\nNOTE: Header cannot contain \"-\" except before the read depth, e.g.:\n";
	print "\t>Name-3
	AGTAGTAGTAG....
Where \"Name\" is the sequence name and \"3\" is the depth.\n\n";
	print "Options:\n";
	print "\t-i	: Path to input file (fasta)\n";
	print "\t-m	: Maximum stack depth to print [default: not set]\n";
	print "\t-n	: Skip clusters with less than \"n\" depth [default: not set]\n";
	print "\t-x	: Skip clusters with more than \"x\" depth [default: not set]\n";
	print "\t-o	: Output file name. [Default = out.phy]\n";
	print "\n\n";
}

#parse arguments
sub parseArgs{

  my( $params ) =  @_;
  my %opts = %$params;
  
  #defaults
  my $in = $opts{i} or die "\nNo input was provided.\n\n";
  my $min = $opts{n} || 0; 
  my $max = $opts{m} || 0; 
  my $cap = $opts{x} || 0;
  my $out = $opts{o} || "out.fasta"; 
  #return
  return ($in, $out, $min, $max, $cap);
}
