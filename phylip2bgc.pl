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
getopts( 'p:i:1:2:a:o:h', \%opts );

# kill if help option is true
if( $opts{h} ){
  &help;
  die "Help menu\n\n";
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
open(ADMIX, "> admix.csv");
open(P1DATA, "> p1data.csv");
open(P2DATA, "> p2data.csv");

#format and print files for bgc
for (my $loc = 0; $loc < $nchar; $loc++){
	
	#write p1 file
	my $pop = 0;
	print P1DATA "locus ",$loc+1, "\n";
	foreach my $id (@pop1){
		print P1DATA "pop ", $pop, "\n";
		foreach my $ind (keys %{$pop1Ref}){
			if ($assignRef->{$ind} eq $id){
				my $nuc = ${$pop1Ref->{$ind}}->[$loc];	
				$nuc =~ s/A/1	1/gi;
				$nuc =~ s/T/2	2/gi;
				$nuc =~ s/G/3	3/gi;
				$nuc =~ s/C/4	4/gi;
				$nuc =~ s/W/1	2/gi;
				$nuc =~ s/R/1	3/gi;
				$nuc =~ s/M/1	4/gi;
				$nuc =~ s/K/3	2/gi;
				$nuc =~ s/Y/2	4/gi;
				$nuc =~ s/S/3	4/gi;
				$nuc =~ s/-/-9	-9/gi;
				$nuc =~ s/N/-9	-9/gi;
				print P1DATA $nuc, "\n";
			}	
		}
		$pop++;
	}


	#Write p2 file
	$pop = 0;
	print P2DATA "locus ",$loc+1, "\n";
	foreach my $id (@pop2){
		print P2DATA "pop ", $pop, "\n";
		foreach my $ind (keys %{$pop2Ref}){
			if ($assignRef->{$ind} eq $id){
				my $nuc = ${$pop2Ref->{$ind}}->[$loc];	
				$nuc =~ s/A/1	1/gi;
				$nuc =~ s/T/2	2/gi;
				$nuc =~ s/G/3	3/gi;
				$nuc =~ s/C/4	4/gi;
				$nuc =~ s/W/1	2/gi;
				$nuc =~ s/R/1	3/gi;
				$nuc =~ s/M/1	4/gi;
				$nuc =~ s/K/3	2/gi;
				$nuc =~ s/Y/2	4/gi;
				$nuc =~ s/S/3	4/gi;
				$nuc =~ s/-/-9	-9/gi;
				$nuc =~ s/N/-9	-9/gi;
				print P2DATA $nuc, "\n";
			}	
		}
		$pop++;
	}


	#Populate admix file
	$pop = 0;
	print ADMIX "locus ",$loc+1, "\n";
	foreach my $id (@popA){
		print ADMIX "pop ", $pop, "\n";
		foreach my $ind (keys %{$popaRef}){
			if ($assignRef->{$ind} eq $id){
				my $nuc = ${$popaRef->{$ind}}->[$loc];	
				$nuc =~ s/A/1	1/gi;
				$nuc =~ s/T/2	2/gi;
				$nuc =~ s/G/3	3/gi;
				$nuc =~ s/C/4	4/gi;
				$nuc =~ s/W/1	2/gi;
				$nuc =~ s/R/1	3/gi;
				$nuc =~ s/M/1	4/gi;
				$nuc =~ s/K/3	2/gi;
				$nuc =~ s/Y/2	4/gi;
				$nuc =~ s/S/3	4/gi;
				$nuc =~ s/-/-9	-9/gi;
				$nuc =~ s/N/-9	-9/gi;
				print ADMIX $nuc, "\n";
			}	
		}
		$pop++;
	}	
	
}

print "Done writing P1DATA file <p1data.csv>\n";
close P1DATA;
print "Done writing P2DATA file <p2data.csv>\n";
close P2DATA;
print "Done writing ADMIX file <admix.csv>\n\n";
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
	print "\t-n	: Prop. of missing data (including gaps) allowed in parent pop to retain SNP [Def = 0.0] \n";
	print "\t-o	: Output file name. [Default = out.phy]\n";
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
        print "Arr: ",join(",",@arr),"\n";
        print "Ref: ",\@arr,"\n";
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
