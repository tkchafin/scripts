#! /usr/bin/perl

# Contributions by Tyler K. Chafin, Steven M. Mussmann, Max R. Bangs
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
getopts( 'p:i:1:2:a:o:hcn:', \%opts );

# kill if help option is true
if( $opts{h} ){
  &help;
  die "Help menu\n\n";
}

my $comb = 0;
if( $opts{c} ){
  $comb = 1;
}
 
#get options 
my ($map, $phy, $p1, $p2, $a, $out) = &parseArgs(\%opts); 

#Extract pops into an array
my @pop1 = split(/\+/,$p1);
my @pop2 = split(/\+/,$p2);
my @popA = split(/\+/,$a);

$p1 = join(', ',@pop1);
$p2 = join(', ',@pop2);
$a = join(', ', @popA); 

#parse popmap file into a hash with ind as key and popID as value
my $assignRef = &parsePopmap($map); 

#parse phylip file into a hash with ind as key and array of seqs as value
my ($allRef, $ntax, $nchar) = &parsePhylip($phy);

#Print argument report
print "\nPopmap file is: $map\n";
print "Phylip file is: $phy\n";
print "Population 1 is: $p1\n";
print "Population 2 is: $p2\n";
print "Admixed population is: $a\n";
print "Total taxa in phylip file: $ntax\n"; 
print "Total characters in phylip data matrix: $nchar\n";

#Get pop alignments only with ind as key and array of seqs as value
my $pop1Ref = &getMultPops($assignRef, $allRef, \@pop1);
my $pop2Ref = &getMultPops($assignRef, $allRef, \@pop2);
my $popaRef = &getMultPops($assignRef, $allRef, \@popA);

#Check if pops contain data
my $num1 = keys %{$pop1Ref};
my $num2 = keys %{$pop2Ref};
my $num3 = keys %{$popaRef};

if ($num1< 1){
  die "Error: No individuals for pop ID <$p1> were found!\n\n";
}elsif ($num2 < 1){
  die "Error: No individuals for pop ID <$p2> were found!\n\n";
}elsif ($num3 < 1){
  die "Error: No individuals for pop ID <$a> were found!\n\n";
}else{
  print "Found <$num1> individuals in population 1\n";
  print "Found <$num2> individuals in population 2\n";
  print "Found <$num3> individuals in admixed population\n\n";
}


#Open filestreams
open(ADMIX, "> bgc_admix.txt");
open(P1DATA, "> bgc_p1in.txt");
open(P2DATA, "> bgc_p2in.txt");

#for each locus
for (my $loc = 0; $loc < $nchar; $loc++){
	
	#get p1 counts
	
	my %p1counts = (
		"A" => 0,
		"T" => 0,
		"G" => 0,
		"C" => 0
	);
	foreach my $id (@pop1){
		foreach my $ind (keys %{$pop1Ref}){
			if ($assignRef->{$ind} eq $id){
				my $nuc = ${$pop1Ref->{$ind}}->[$loc];	
				if ($nuc eq "A") { $p1counts{"A"} += 2;} #allele count ONLY counts acceptable diploid genotypes
				elsif ($nuc eq "T") { $p1counts{"T"} += 2;}
				elsif ($nuc eq  "C") { $p1counts{"C"} += 2;}
				elsif ($nuc eq  "G") { $p1counts{"G"} += 2;}
				elsif ($nuc eq  "R") { $p1counts{"A"}++; $p1counts{"G"}++;}
				elsif ($nuc eq  "Y") { $p1counts{"C"}++; $p1counts{"T"}++;}
				elsif ($nuc eq  "S") { $p1counts{"G"}++; $p1counts{"C"}++;}
				elsif ($nuc eq  "W") { $p1counts{"A"}++; $p1counts{"T"}++;}
				elsif ($nuc eq  "K") { $p1counts{"G"}++; $p1counts{"T"}++;}
				elsif ($nuc eq  "M") { $p1counts{"C"}++; $p1counts{"A"}++;}
			}	
		}
	}

	#get p2 counts
	my %p2counts = (
		"A" => 0,
		"T" => 0,
		"G" => 0,
		"C" => 0
	);
	foreach my $id (@pop2){
		foreach my $ind (keys %{$pop2Ref}){
			if ($assignRef->{$ind} eq $id){
				my $nuc = ${$pop2Ref->{$ind}}->[$loc];	
				if ($nuc eq "A") { $p2counts{"A"} += 2;}
				elsif ($nuc eq "T") { $p2counts{"T"} += 2;}
				elsif ($nuc eq  "C") { $p2counts{"C"} += 2;}
				elsif ($nuc eq  "G") { $p2counts{"G"} += 2;}
				elsif ($nuc eq  "R") { $p2counts{"A"}++; $p2counts{"G"}++;}
				elsif ($nuc eq  "Y") { $p2counts{"C"}++; $p2counts{"T"}++;}
				elsif ($nuc eq  "S") { $p2counts{"G"}++; $p2counts{"C"}++;}
				elsif ($nuc eq  "W") { $p2counts{"A"}++; $p2counts{"T"}++;}
				elsif ($nuc eq  "K") { $p2counts{"G"}++; $p2counts{"T"}++;}
				elsif ($nuc eq "M") { $p2counts{"C"}++; $p2counts{"A"}++;}
			}	
		}
	}

	#choose first and second allele
	#NOTE: Script will SKIP loci that are monomorphic in parental pops, or >2 alleles
	my $first = "";
	my $second = "";
	for my $k1 (keys %p1counts){
		if ($p1counts{$k1} > 0){
			if ($first eq ""){
				$first = $k1;
			}elsif($second eq ""){
				$second = $k1;
			}else{
				#print $first, "\t", $second, "---";
				print("Locus $loc not bi-allelic. Skipping.\n");
				next;
			}
		}
	}
	# if ($first eq $second){
	# 	print $first, "\t", $second, "---";
	# 	print "Something went wrong...\n";
	# 	exit;
	# }
	
	if ($first eq ""){
		print("Locus $loc had no alleles in parental population 1. Skipping.\n");
		next;
	}
	for my $k2 (keys %p2counts){
		if ($p2counts{$k2} > 0){
			if($first ne $k2 && $second eq ""){
				$second = $k2;
			}elsif($second eq $k2 || $first eq $k2){
				next;
			}else{
				#print $first, "\t", $second, " --- ", $k2," --- ";
				print("Locus $loc not bi-allelic. Skipping.\n");
				next;
			}
		}
	}
	if ($second eq ""){
		print("Locus $loc monomorphic in parental populations. Skipping.\n");
		next;
	}

	#write to parental files
	print P1DATA "locus ",$loc, "\n";
	print P1DATA $p1counts{$first}, "\t", $p1counts{$second}, "\n";
	print P2DATA "locus ",$loc, "\n";
	print P2DATA $p2counts{$first}, "\t", $p2counts{$second}, "\n";
	

	#Populate admix file
	print ADMIX "locus ",$loc, "\n";
	if ($comb == 0){
		my $pop = 0;
		foreach my $id (@popA){
			print ADMIX "pop ", $pop, "\n";
			foreach my $ind (keys %{$popaRef}){
				if ($assignRef->{$ind} eq $id){
					my %indcounts = (
						"A" => 0,
						"T" => 0,
						"G" => 0,
						"C" => 0
					);
					my $nuc = ${$popaRef->{$ind}}->[$loc];	
					if ($nuc eq "A") { $indcounts{"A"} += 2;}
					elsif ($nuc eq "T") { $indcounts{"T"} += 2;}
					elsif ($nuc eq  "C") { $indcounts{"C"} += 2;}
					elsif ($nuc eq  "G") { $indcounts{"G"} += 2;}
					elsif ($nuc eq  "R") { $indcounts{"A"}++; $indcounts{"G"}++;}
					elsif ($nuc eq  "Y") { $indcounts{"C"}++; $indcounts{"T"}++;}
					elsif ($nuc eq  "S") { $indcounts{"G"}++; $indcounts{"C"}++;}
					elsif ($nuc eq  "W") { $indcounts{"A"}++; $indcounts{"T"}++;}
					elsif ($nuc eq  "K") { $indcounts{"G"}++; $indcounts{"T"}++;}
					elsif ($nuc eq "M") { $indcounts{"C"}++; $indcounts{"A"}++;}
					
					#write genotype for this individual
					if(($indcounts{$first} + $indcounts{$second}) != 2){
						print ADMIX "-9\t-9\n";
					}else{
						print ADMIX $indcounts{$first}, "\t", $indcounts{$second}, "\n";
					}
				}	
			}
			$pop++;
		}
	}else{
		print ADMIX "locus ",$loc, "\n";
		print ADMIX "pop 0\n";
		foreach my $id (@popA){
			foreach my $ind (keys %{$popaRef}){
				if ($assignRef->{$ind} eq $id){
					my %indcounts = (
						"A" => 0,
						"T" => 0,
						"G" => 0,
						"C" => 0
					);
					my $nuc = ${$popaRef->{$ind}}->[$loc];	
					if ($nuc eq "A") { $indcounts{"A"} += 2;}
					elsif ($nuc eq "T") { $indcounts{"T"} += 2;}
					elsif ($nuc eq  "C") { $indcounts{"C"} += 2;}
					elsif ($nuc eq  "G") { $indcounts{"G"} += 2;}
					elsif ($nuc eq  "R") { $indcounts{"A"}++; $indcounts{"G"}++;}
					elsif ($nuc eq  "Y") { $indcounts{"C"}++; $indcounts{"T"}++;}
					elsif ($nuc eq  "S") { $indcounts{"G"}++; $indcounts{"C"}++;}
					elsif ($nuc eq  "W") { $indcounts{"A"}++; $indcounts{"T"}++;}
					elsif ($nuc eq  "K") { $indcounts{"G"}++; $indcounts{"T"}++;}
					elsif ($nuc eq "M") { $indcounts{"C"}++; $indcounts{"A"}++;}
					
					#write genotype for this individual
					if(($indcounts{$first} + $indcounts{$second}) != 2){
						print ADMIX "-9\t-9\n";
					}else{
						print ADMIX $indcounts{$first}, "\t", $indcounts{$second}, "\n";
					}
				}	
			}
		}
	}
}

print "Done writing P1DATA file <bgc_p1in.txt>\n";
close P1DATA;
print "Done writing P2DATA file <bgc_p2in.txt>\n";
close P2DATA;
print "Done writing ADMIX file <bgc_admix.txt>\n\n";
close ADMIX;


 ########################### SUBROUTINES ###############################

 sub help{
	 
	print "\nphylip2introgress.pl by Tyler Chafin with some code from Max Bangs and Steve Mussmann\n";
	print "\nThis script converts from phylip format to the input for bgc.\n";
	print "A population map should be given in a tab-delimited file, formatted as:\n";
	print "\n\tSampleName\tPopID\n\n";
	print "Where PopID can be a string or integer, and SampleName must exactly match a corresponsing line in the phylip file. Any samples not in the popmap will not be included in the output files.\n\n";
	print "Options:\n";
	print "\t-p	: Path to popmap file (tab-delimited)\n";
	print "\t-1	: Identifier for population 1 (include multiple as: pop1+pop2)\n";
	print "\t-2	: Identifier for population 2 (include multiple as: pop1+pop2)\n";
	print "\t-a	: Identifier for admixed population(s) (include multiple as: pop1+pop2)\n";
	print "\t-i	: Path to input file (phylip)\n";
	print "\t-c	: Combine admixed populations under a single label [boolean]\n";
	print "\t-h	: Displays this help message\n";
	print "\n\n";
}

#parse arguments
sub parseArgs{

  my( $params ) =  @_;
  my %opts = %$params;
  
  #defaults
  my $map = $opts{p} or die "\nPopmap not specified.\n\n";
  my $phy = $opts{i} or die "\nPhylip file not specified.\n\n";
  my $p1  = $opts{1} or die "\nPopulation 1 not specified.\n\n";
  my $p2  = $opts{2} or die "\nPopulation 2 not specified.\n\n";
  my $a  = $opts{a} or die "\nNo admixed population specified.\n\n";
  my $out = $opts{o} || "out.phy"; 
  #return
  return ($map, $phy, $p1, $p2, $a, $out);
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

#parse phylip file -> This version returns array refs, not strings, of sequences
sub parsePhylip{
	
	my $toParse = shift(@_); 
	
	#vars
	my %toReturn; 
	my @seq; 
	my $ntax;
	my $nchar;
	
	#open popmap
	open (PHY, $toParse) or die "Oh no! Cannot open $toParse: $!\n"; 
	
	my $num = 0; 
	while (my $line = <PHY>){
	  $num++; 
	  chomp $line; 
	  
	  #Skip first line
	  if ($num == 1){
		my @temp = split( /\s+/, $line);
		$ntax = $temp[0];
		$nchar = $temp[1];
	    next; 
	  }
	  
	  #ignore if blank
      if( $line =~ /^\w/ ){
        my @temp = split( /\s+/, $line); 
        my @arr = split(//, $temp[1]);
        #push array ref into our hash
        $toReturn{$temp[0]} = \@arr;
      }
	}
	
	close PHY;
	return( \%toReturn, $ntax, $nchar);
}

#Get alignments for only populations of interest 
sub getPops{
	my $pops = $_[0];
	my $seqs = $_[1];
	my $first = $_[2];
	my $second = $_[3];
	
	my %pop1;
	my %pop2;
	
	foreach my $key (keys %{$pops}){
		#If pop ID matches
		if ($pops->{$key} eq $first){
			${$pop1{$key}} = $seqs->{$key}; 
		}elsif ($pops->{$key} eq $second){
			${$pop2{$key}} = $seqs->{$key}; 
		}
	}
	return(\%pop1, \%pop2);
}

#Modified getPops subroutine, gets all pops matching array of options, returns as one hash
sub getMultPops{
	
	my $popRef 		= $_[0];
	my $seqRef 		= $_[1];
	my $toGetRef 	= $_[2];
	
	my %pop; 
	
	foreach my $id (@{$toGetRef}){
		foreach my $key (keys %{$popRef}){
			#If pop ID matches, get sequence
			if ($popRef->{$key} eq $id){
				if (exists $seqRef->{$key}){
					${$pop{$key}} = $seqRef->{$key}; 
				}else{
					delete $popRef->{$key}; #Remove entry from pop hash
					print "Warning: Sample <$key> was not found in sequence file. Skipping.\n";
				}
			}
		}
	}
	return(\%pop);	
}
	

 
# subroutine to print data out to a phylip file

sub phyprint{

  my( $out, $hashref ) = @_;
  
  # get the number of sequences
  my $seqs = scalar keys %$hashref;

  # get the length of the alignment
  my $alignlength;
  foreach my $key( sort keys %{ $hashref } ){
    $alignlength = length( ${$hashref}{$key} );
  }
  
  # get the length of the longest 
  my $keylength = 0;
  foreach my $key( sort keys %{ $hashref } ){
    my $temp = length( $key );
    if( $temp > $keylength ){
      $keylength = $temp;
    }
  }

  # open the output file for printing
  open( OUT, '>', $out ) or die "Can't open $out, d-bag: $!\n\n";

  # print the first line to the phylip file
  print OUT "$seqs $alignlength\n";

  # print the hash
  foreach my $key( sort keys %{ $hashref } ){
    my $templength = length( $key );
    my $whitespace = ( ( $keylength + 2 ) - $templength );
    print OUT $key;
    for( my $i=0; $i<$whitespace; $i++ ){
      print OUT " ";
    }
    print OUT ${$hashref}{$key}, "\n";
  }

  # close the output file
  close OUT;

}
