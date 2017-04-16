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
getopts( 'g:f:s:x:y:o:ht:r:1:2:3:4:5:6:7:8:9:A:B:C:b', \%opts );

# kill if help option is true
if( $opts{h} ){
  &help;
  die "Help menu\n\n";
}
 
#get options 
my ($geo, $seq, $samp, $x, $y, $out, $skip, $thresh, $rand, $resample) = &parseArgs(\%opts); 
my ($A, $B, $C, $a1, $a2, $a3, $a4, $a5, $a6, $a7, $a8, $a9) = &parseSAMOVAArgs(\%opts); 

#Get FASTA sequences, read into hash
my %seqHash; 
my $sample; 
my $getNucsRef = ""; 
open (ALN, "$seq") || die "Error: Cannot open file $seq: $!\n\n"; 
while (<ALN>){
  chomp; 
  if ($_ =~ /^>/){
    $sample = $_; 
    $sample =~ s/>//;
    next;
  }else{
    if ($rand != 0){
      if ($getNucsRef eq ""){ #If first sample, generate list of positions to sample
		open (RAN, ">$out.sampled") || die "Error: Cannot open file $out.sampled: $!\n";
		if ($rand > length($_)){
			my $len = length($_)+1;
			print "Warning! You asked to sample more columns than exist in your alignment! Setting <-r> to sequence length ($len)\n";
		}
        $getNucsRef = &generatePositions(length($_), $rand, $resample); 
      } 
      #Randomly sample sequence columns
	  my $newString = &sampleString($_, $getNucsRef);
	  print RAN ">$sample\n$newString\n";
	  $seqHash{$sample} = $newString;
      $sample = "";
      next;
    }else{
      $seqHash{$sample} = $_;
      #print "$sample: $_\n";
      $sample = "";
      next;
    }
  }
}
close ALN; 
$rand and close RAN;


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
if ($rand){
  print "Number of nucleotide columns randomly sampled: $rand\n";
  my $bool = "FALSE";
  $resample and $bool = "TRUE";
  print "Nucleotides sampled with replacement: $bool\n";
  print "Outputted resampled alignment as FASTA to: $out.sampled\n";
}
print "\nDone!\n\n";
exit;

 ########################### SUBROUTINES ###############################

 sub help{
	 
	print "\nmakeSAMOVA.pl by Tyler K. Chafin; tkchafin\@uark.edu\n";
	print "\nThis script takes FASTA sequence alignment and creates the input for SAMOVA 2.0\n\n";
	print "A tab-delimited file containing columns for sample names and geographic coordinates should be provided. By default, sample names are assumed to be in column 1, with X and Y coordinates in columns 2 and 3.\n";
	print "\nNOTE: Sample names should correspond exactly across sequence and coordinates files.\n\n";
	print "NOTE: Sample clustering is in a very primitive state. A hard numerical 'distance' will be calculated as the hypotenuse length using X and Y distances (calculated as the absolute value of difference in coordinate values), and if below the threshold value will be clumped. Population coordinate will be determined using the first sampled individual for the population.\n\n";
	print "NOTE: This script will output with WINDOWS EOL markers\n\n";
	print "Mandatory arguments:\n";
	print "\t-g	: Path to coordinates file (tab-delimited)\n";
	print "\t-f	: Path to sequence file (fasta)\n\n";
	print "Optional arguments:\n";
	print "\t-t	: Threshold difference in coordinates to cluster samples into a pop [default=0.001]\n";
	print "\t-s	: Column in coordinates file with sample IDs [default=1]\n";
	print "\t-x	: Column in coordinates file with X coordinate [default=2]\n";
	print "\t-y	: Column in coordinates file with Y coordinate [default=3]\n";
	print "\t-o	: Output file prefix [Default = out]\n";
	print "\t-k : Lines to skip in coordinates file [default = 1, assumes header line]\n";
	print "\t-r : Randomly sample X positions from sequences [Default=0, no sampling]\n";
	print "\t-b : Toggle on replacement for random sampling- can be used to generate bootstraps\n";
	print "\t-h	: Displays this help message\n";
	
	print "\nSAMOVA arguments:\n";
	print "\t-1	: DISTMETAMOVA [default=5]\n";
	print "\t-2	: GAMMAVALUE [default=0.17]\n";
	print "\t-3	: MISSING_LEVEL [default=0.05]\n";
	print "\t-4	: TRANSITION_WEIGHT [default=1]\n";
	print "\t-5	: TRANSVERSION_WEIGHT [default=1]\n";
	print "\t-6	: DELETION_WEIGHT [default=1]\n";
	print "\t-7	: NUM_GROUPS [default=2]\n";
	print "\t-8	: NUM_INIT_CONDS [default=10]\n";
	print "\t-9	: NO_GEO_STRUCTURE [default=0]\n";
	print "\t-A	: SA_NUM_STEPS [default=10000]\n";
	print "\t-B	: SA_EXP_FACTOR [default=0.915811421]\n";
	print "\t-C	: SIMULATIONS [default=0]\n";
	
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
  my $r = $opts{r} || 0;
  my $b = $opts{b} || 0;
  #return
  return ($geot, $seqt, $st, $xt, $yt, $outt, $skipt, $threshT, $r, $b);
}

#parse SAMOVA arguments
sub parseSAMOVAArgs{

  my( $params ) =  @_;
  my %opts = %$params;
  
  #defaults
  my $t1 = $opts{1} || 5;
  my $t2 = $opts{2} || 0.17;
  my $t3 = $opts{3} || 0.05;
  my $t4 = $opts{4} || 1;
  my $t5 = $opts{5} || 1;
  my $t6 = $opts{6} || 1;
  my $t7 = $opts{7} || 2;
  my $t8 = $opts{8} || 10;
  my $t9 = $opts{9} || 0;
  my $ta = $opts{A} || 10000;
  my $tb = $opts{B} || 0.915811421;
  my $tc = $opts{B} || 0;
  #return
  return ($ta, $tb, $tc, $t1, $t2, $t3, $t4, $t5, $t6, $t7, $t8, $t9);
}


#generate list of integer positions to sample
sub generatePositions{
	#print "From generatePositions:\n";
	my $max = $_[0];
	my $howMany = $_[1];
	my $how = $_[2];
	my %returnMap; 
	my @returnArray;
	for (my $i=0; $i < $howMany; $i++){
		#$how == 0 and my $num = int(&randomNumsWithoutReplacement($max, \%returnMap));
		my $num = int(&randomNumsWithoutReplacement($max, \%returnMap));
		#$how == 1 and my $num = int(&randomNumsWitReplacement($max));
		#print "Found a number: $num\n";
		$returnMap{$num} = "";
		push(@returnArray, $num);
	}
	return \@returnArray;
}

sub randomNumsWithoutReplacement{
	my $m = $_[0];
	my $already = $_[1];
	my $returnNum = rand($m); 
	while(exists $already->{$returnNum}){
		$returnNum = rand($m); 
	}
	return $returnNum;
}

sub randomNumsWithReplacement{
	my $m = $_[0];
	my $returnNum = rand($m); 
	return $returnNum;
}

sub sampleString{

	my $oldString = $_[0];
	my $toSample = $_[1];
	my $returnString = ""; 

	for my $i(0 .. $#{$toSample}){
		$returnString .= substr($oldString, ($toSample->[$i])-1, 1);
	}
	return $returnString;
}








