#! /usr/bin/perl

# Written by: Tyler K. Chafin
# Contact: tkchafin@uark.edu

use strict; 
use warnings; 
use Getopt::Std; 

#Die if no arguments given
if( scalar( @ARGV ) == 0 ){
  &help;
  die "No options given\n\n";
}

#Parse arguments
my %opts;
getopts( 'g:f:s:x:y:o:ht:', \%opts );

# kill if help option is true
if( $opts{h} ){
  &help;
  die "Help menu\n\n";
}
 
#get options 
my ($geo, $seq, $samp, $x, $y, $out, $skip, $thresh) = &parseArgs(\%opts); 

#Get FASTA sequences, read into hash
my %seqHash; 
my $sample; 
open (ALN, "$seq") || die "Error: Cannot open file $seq: $!\n\n"; 
while (<ALN>){
  chomp; 
  if ($_ =~ /^>/){
    $sample = $_; 
    $sample =~ s/>//;
    next;
  }else{
    $seqHash{$sample} = $_;
    #print "$sample: $_\n";
    $sample = "";
    next;
  }
}
close ALN; 

#Get coordinates from coordinate file for each sample
open (GEO, "$geo") || die "Error: Cannot open file $geo: $!\n\n";
my %geoHash;
my $count = 0;
while (<GEO>){
  $count++;
  $count <= $skip and next;
  chomp; 
  my @temp = split /\t/, $_;
  my $name = $temp[$samp-1];
 # print "X is $temp[$x-1] and Y is $temp[$y-1]\n";
  my $coordString = $temp[$x-1];
  $coordString .= "\t" . $temp[$y-1];
  $geoHash{$name} = $coordString;
}
close GEO;

#Loop through and cluster populations by identical coordinates
#Extend in future to clump within a threshold? 
my %popHash; 
foreach my $key (sort keys %geoHash){
  #print "$key is $geoHash{$key}\n";
  #If popHash is empty, add sample 
  #print scalar(keys(%popHash)) . "\n";
  if (scalar(keys %popHash) <= 0){
    my @temp; 
    print "Population hash empty. Initializing population for sample $key\n";
    $temp[0] = $geoHash{$key}; 
    $temp[1] = $key; 
    $popHash{$key} = [@temp];
    next;
  }
  my $added = 0;
  foreach my $pop (sort keys %popHash){
    #if sample has near coords to existing pop, add
    #otherwise, make new pop
    #print "Check if $geoHash{$key} matches $popHash{$pop}->[0]"; 

    if ($geoHash{$key} eq $popHash{$pop}->[0]){
	  #print "..same\n";
      push(@{$popHash{$pop}}, $key);
      $added=1; 
      print "Match:\tDistance between $key and $pop: is 0.000!\n";
      last;
    }else{
      my @check = split /\t/, $geoHash{$key};
      my @compare = split /\t/, $popHash{$pop}->[0]; 
      my $xdist = abs($check[0] - $compare[0]);
      my $ydist = abs($check[1] - $compare[1]);
      my $dist = abs(sqrt(($xdist**2)+($ydist**2)));
      #print "Distance between $key and $pop: is $dist!\n";
      if ($dist <= $thresh){
		#print "Distance $dist is less than threshold $thresh\n";
		print "Match:\tDistance between $key and $pop: is $dist!\n";
        push(@{$popHash{$pop}}, $key);
        $added=1; 
        last;
	  }
    }  
  }	
  if ($added == 0){ #Sample not found anywhere
  my @temp; 
    print "Sample $key did not match. Starting new population.\n";
    $temp[0] = $geoHash{$key}; 
    $temp[1] = $key; 
    $popHash{$key} = [@temp];
  }
}
#print scalar(keys(%popHash)) . "\n";

#Print output files
open (COORD, ">$out.geo") || die "Error: Cannot open file $out.geo: $!\n";
open (ARP, ">$out.arp") || die "Error: Cannot open file $out.arp: $!\n";

print ARP "#ARlequin input file written by makeSAMOVA.pl\n\n";
print ARP "[PROFILE]\n\tTitle=\"SAMOVA\"\n";
my $nb = keys %popHash;
print ARP "\tNbSamples=$nb\n\n\tGenotypicData=0\n\tGameticPhase=0\n";
print ARP "\tRecessiveData=0\n\tDataType=DNA\n\tLocusSeparator=NONE\n";
print ARP "\tMissingData='?'\n\n[Data]\n\t[[Samples]]\n\n";

my $popCount=1;
my $singletonCount=0;
foreach my $pop (sort keys %popHash){
  print COORD "$popCount\t\"pop$pop\"\t$popHash{$pop}->[0]\t1\n";
  print ARP "\t\tSampleName=\"pop$pop\"\n";
  my $size = $#{$popHash{$pop}};
  print ARP "\t\tSampleSize=$size\n\t\tSampleData={\n";
  $size == 1 and $singletonCount++;
  for (my $i=1; $i <= $#{$popHash{$pop}}; $i++){
    print ARP "$popHash{$pop}->[$i]\t1\t$seqHash{$popHash{$pop}->[$i]}\n";
  }
  print ARP "\n}\n";
  $popCount++;
}
close COORD;
close ARP;

#Report some stuff
print "\nArlequin input written to: $out.arp\n";
print "Coordinate outputs written to: $out.geo\n";
print "Number of collapsed populations: $popCount\n";
print "Number of populations of only one individual: $singletonCount\n";
$singletonCount>0 and print "Warning: SAMOVA may not work properly with singleton sites!\n\n";
print "Done!\n\n";
exit;

 ########################### SUBROUTINES ###############################

 sub help{
	 
	print "\nmakeSAMOVA.pl by Tyler K. Chafin; tkchafin\@uark.edu\n";
	print "\nThis script takes FASTA sequence alignment and creates the input for SAMOVA 2.0\n\n";
	print "A tab-delimited file containing columns for sample names and geographic coordinates should be provided. By default, sample names are assumed to be in column 1, with X and Y coordinates in columns 2 and 3.\n";
	print "\nNOTE: Sample names should correspond exactly across sequence and coordinates files.\n\n";
	print "NOTE: Sample clustering is in a very primitive state. A hard numerical 'distance' will be calculated as the hypotenuse length using X and Y distances (calculated as the absolute value of difference in coordinate values), and if below the threshold value will be clumped. Population coordinate will be determined using the first sampled individual for the population.\n\n";
	
	print "Mandatory arguments:\n";
	print "\t-g	: Path to coordinates file (tab-delimited)\n";
	print "\t-f	: Path to sequence file (fasta)\n\n";
	print "Optional arguments:\n";
	print "\t-s	: Threshold difference in coordinates to cluster samples into a pop [default=0.001]\n";
	print "\t-s	: Column in coordinates file with sample IDs [default=1]\n";
	print "\t-x	: Column in coordinates file with X coordinate [default=2]\n";
	print "\t-y	: Column in coordinates file with Y coordinate [default=3]\n";
	print "\t-o	: Output file prefix [Default = out]\n";
	print "\t-k : Lines to skip in coordinates file [default = 1, assumes header line]\n";
	print "\t-h	: Displays this help message\n";
	print "\n\n";
}

#parse arguments
sub parseArgs{

  my( $params ) =  @_;
  my %opts = %$params;
  
  #defaults
  my $geot = $opts{g} or die "\nCoordinates file not specified.\n\n";
  my $seqt = $opts{f} or die "\nFasta file not specified.\n\n";
  my $st  = $opts{s} || 1;
  my $xt  = $opts{x} || 2;
  my $yt  = $opts{y} || 3;
  my $outt = $opts{o} || "out"; 
  my $skipt = $opts{k} || 1;
  my $threshT = $opts{t} || 0.001; 
  #return
  return ($geot, $seqt, $st, $xt, $yt, $outt, $skipt, $threshT);
}

#parse popmap file
sub parsePopmap{
	
	my $toParse = $_[0]; 
	
	#vars
	my %toReturn; 
	
	#open popmap
	open (POP, $toParse) or die "Oh no! Cannot open $toParse: $!\n"; 
	
	while (my $line = <POP>){
	  chomp $line; 

	  #ignore if blank
      if( $line =~ /^\w/ ){
        my @temp = split( /\s+/, $line); 

        #push into our hash
        $toReturn{$temp[0]} = $temp[1];
      }
	}
	
	close POP;
	return( \%toReturn);
}
