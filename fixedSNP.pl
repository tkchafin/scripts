#! /usr/bin/perl

# By Tyler K. Chafin
# Contact: tkchafin@uark.edu

use strict; 
use warnings; 
use Getopt::Std; 

#Die if no arguments given
if( scalar( @ARGV ) == 0 ){
  &help;
  die "Y u no give option???\n\n";
}

#Parse arguments
my %opts;
getopts( 'p:i:hn:1:2:g:o:x', \%opts );

# kill if help option is true
if( $opts{h} ){
  &help;
  die "U want halp?\n\n";
}

my $skip = 0;
if ($opts{x}){
  $skip = 1; 
  print "Warning: You have chosen to skip checking for fixed SNPs...\n";
  print "...Only filtering based on N and gap content.\n";
} 
 
#get options 
my ($map, $phy, $n, $p1, $p2, $gap, $out) = &parseArgs(\%opts); 

print "\nPopmap file is: $map\n";
print "Phylip file is: $phy\n";
print "Population 1 is: $p1\n";
print "Population 2 is: $p2\n";
print "Gap threshold is set to: $gap\n";
print "N threshold is set to: $n\n";

#parse popmap file into a hash with ind as key and popID as value
my $assignRef = &parsePopmap($map); 

#parse phylip file into a hash with ind as key and array of seqs as value
my $allRef = &parsePhylip($phy);

#Get pop1 and pop2 alignments only with ind as key and array of seqs as value
my ($pop1Ref, $pop2Ref) = &getPops($assignRef, $allRef, $p1, $p2); 

#Check if pops contain data
my $num1 = keys %{$pop1Ref};
my $num2 = keys %{$pop2Ref};

if ($num1< 1){
  die "Error: No individuals for $p1 were found!\n\n";
}elsif ($num2 < 1){
  die "Error: No individuals for $p2 were found!\n\n";
}else{
  print "Found <$num1> individuals in population 1\n";
  print "Found <$num2> individuals in population 2\n";
}

#filter for N content
#parasitized subroutine steve wrote
my $pop1Align = &getColumns($pop1Ref);
my $pop2Align = &getColumns($pop2Ref);

#Parse pop1 and pop2 alignments for fixed SNPs
my $toDelete = &parsePopAlignment($pop1Align, $pop2Align, $n, $gap, $skip);

#Delete flagged columns from the full alignment
#another one parasitized from steve
removeColumns($allRef, $toDelete);

#print final result, also subroutine parasitized from Steve
phyprint($out, $allRef);

#Report matrix content
print "\n";
reportMatrixContent($allRef);

 ########################### SUBROUTINES ###############################

 sub help{
	 
	print "\nfixedSNP.pl is a perl script developed by Tyler K. Chafin\n";
	print "\nThis script filters a phylip alignment to only include SNPs which are";
	print " differentially fixed for two given populations. Population assignments";
	print " should be given in a tab-delimited file, formatted as:\n";
	print "\n\tSampleName\tPopID\n\n";
	print "Where PopID can be a string or integer, and SampleName must exactly match";
	print " a corresponsing line in the phylip file\n\n";
	print "Options:\n";
	print "\t-p	: Path to popmap file (tab-delimited)\n";
	print "\t-1	: Identifier for population 1\n";
	print "\t-2	: Identifier for population 2\n";
	print "\t-i	: Path to input file (phylip)\n";
	print "\t-n	: Proportion of N's allowed within population to retain a column 
			[Default = 0.5] \n";
	print "\t-g	: Proportion of gaps allowed within population to retain a column 
			[Default = 0.0; e.g. no gaps allowed] \n";
	print "\t-o	: Output file name. [Default = out.phy]\n";
	print "\t-h	: Displays this help message\n";
	print "\t-x : Toggle to skip check for fixed SNP (and only filter Ns)\n";
	print "\n\n";
}

#parse arguments
sub parseArgs{

  my( $params ) =  @_;
  my %opts = %$params;
  
  #defaults
  my $map = $opts{p} or die "\nWhere your popmap at!??\n\n";
  my $phy = $opts{i} or die "\nY u no give Phylip file? You want find SNP or not?!?\n\n";
  my $n   = $opts{n} || 0.5; 
  my $p1  = $opts{1} or die "\nU don't give population 1 wat u thinking\n\n";
  my $p2  = $opts{2} or die "\nwhere ur population 2 at!??\n\n";
  my $gap = $opts{g} || 0.0;
  my $out = $opts{o} || "out.phy"; 
  #return
  return ($map, $phy, $n, $p1, $p2, $gap, $out);
}

#parse popmap file
sub parsePopmap{
	
	my $toParse = $_[0]; 
	
	#vars
	my %toReturn; 
	
	#open popmap
	open (POP, $toParse) or die "Oh NO! Cannot open $toParse: $!\n"; 
	
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

#parse popmap file
sub parsePhylip{
	
	my $toParse = shift(@_); 
	
	#vars
	my %toReturn; 
	my @seq; 
	
	#open popmap
	open (PHY, $toParse) or die "Oh NO! Can no open $toParse: $!\n"; 
	
	my $num = 0; 
	while (my $line = <PHY>){
	  $num++; 
	  chomp $line; 
	  
	  #Skip first line
	  if ($num == 1){
	    next; 
	  }
	  
	  #ignore if blank
      if( $line =~ /^\w/ ){
        my @temp = split( /\s+/, $line); 
        
        #push array ref into our hash
        $toReturn{$temp[0]} = $temp[1];
      }
	}
	
	close PHY;
	return( \%toReturn);
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
	
# subroutine to put sequence alignment into a hash with the index value of the alignment as the key and a string of nucleotides at that index as the value
# modified from a subroutine steve wrote
sub getColumns{

  my( $hashref ) = @_;
  
  my %align; # hash of arrays to hold position in alignment (hash key) and all characters at that position in alignment (value)

  #For each individual
  foreach my $key( sort keys %{ $hashref } ){
    my $index = 0;
    my @seq = split( //, ${$hashref->{$key}}  );
    #for each nucleotide 
    foreach my $item( @seq ){
      $align{$index} .= $item;
      $index++;
    }
  }
  
  return( \%align );

}	
	
#Subroutine to parse the alignment 
sub parsePopAlignment{
	
	my $p1 = $_[0];
	my $p2 = $_[1];
	my $thresholdN = $_[2];
	my $thresholdG = $_[3];
	my $x = $_[4];
	my @blacklist; 
	
	#To track fixed alleles in each pop
	my $alleles1 = parseColumn($p1, $thresholdN, $thresholdG, \@blacklist, $x); 
	my $alleles2 = parseColumn($p2, $thresholdN, $thresholdG, \@blacklist, $x);
	
	#Make sure both pops have same number of columns
	if ((scalar(@{$alleles1})) != (scalar(@{$alleles1}))){
		die "\nError: Y ur populations have not same sequence leNGTH???\n\n";
	}else{
		#Only keep loci which are differentially fixed 
		#Make sure to check anything fixed in pop1 is different 
		#from fixed in pop2 
		for(my $i=0; $i < scalar(@{$alleles1}); $i++){
			my $check1 = $alleles1->[$i] =~ tr/NBDHV-/NBDHV-/;
			my $check2 = $alleles2->[$i] =~ tr/NBDHV-/NBDHV-/;
			#If either pop was variable, or fixed for gaps or Ns
			if ($check1 > 0 || $check2 > 0){
				next; 
			}else{
				#If both fixed for same allele
				if ($alleles1->[$i] eq $alleles2->[$i]){
					push(@blacklist, $i) unless $x; 
					next;
				}
			} 
		} 
	}	
	return(\@blacklist);
}	
	
#internal subroutine to parse columns
sub parseColumn{
	
	my $hash = $_[0];
	my $thresholdN = $_[1];
	my $thresholdG = $_[2];
	my $blistRef = $_[3];
	my $x = $_[4];
	
	my @alleles;
	
	foreach my $key (sort keys %{$hash}){
		#Get length of aligned column
		my $length = length(${$hash}{$key})*2; #Assumes 2N alleles, diplod
		#Get nucleotide frequencies for column
		my ($freqRef) = &parseDipSeq(${$hash}{$key});
		
		#parse frequencies 
		my $ncontent = ($freqRef->{"N"} / $length);
		my $gcontent = ($freqRef->{"-"} / $length);
		
	
		#if gap freq > gap threshold add to blacklist
		if (sprintf( "%.2f", $gcontent ) > sprintf( "%.2f", $thresholdG)){
			#Put N in alleles array and add column to blacklist
			$alleles[$key] = "-";
			push(@{$blistRef}, $key);
			next;
		#if N freq > N threshold add to blacklist	
		}elsif (sprintf( "%.2f", $ncontent ) > sprintf( "%.2f", $thresholdN)){
			#Put N in alleles array and add column to blacklist
			$alleles[$key] = "N";
			push(@{$blistRef}, $key);
			next;
		}else{
			#If any freq = length *2 (assuming diploid)
			my $seen = 0; 
			my $fixed = 0;
			my $pick = "";
			foreach my $nuc(keys %{$freqRef}){
				#Skip Ns and gaps
				if ($nuc eq "N" || $nuc eq "-"){
					next;
				#If nuc has non-zero frequency
				}elsif ($freqRef->{$nuc} != 0){
					if ($seen == 1){
						$fixed = 1;
						last; 
					}else{
						$seen = 1; 
						$pick = $nuc;
						next;
					}
					$pick = $nuc; 
					next; 
				}	
			}
			
			#If a fixed SNP was found
			if ($fixed == 0){		
				#Set allele to V and blacklist column
				$alleles[$key] = $pick;
			}else{
				#SNP must have been variable
				$alleles[$key] = "V";
				push(@{$blistRef}, $key) unless $x; 
				next;
			}
		}
	}
	return(\@alleles);
}
	
#Subroutine to get number of each nucleotide in an alignment
#Assumed diploidy
sub parseDipSeq{
	
	my $seq = shift(@_);
	#Assumes diploidy. iupac codes for 3 nuc ambiguity are treated as Ns (see sub percentn)
	my %codeHash =( 
		"R" => ["A", "G"], 
		"Y" => ["T", "C"], 
		"K" => ["G", "T"], 
		"M" => ["A", "C"], 
		"S" => ["G", "C"],
		"W" => ["A", "T"], 
		"A" => ["A", "A"], 
		"T" => ["T", "T"], 
		"G" => ["G", "G"], 
		"C" => ["C", "C"], 
		"-" => ["-", "-"],
	); 
	#hash to store snp frequencies
	my %temp =(
		"A" => 0,
		"G"	=> 0,
		"C"	=> 0,
		"T"	=> 0,
		"N"	=> 0,
		"-"	=> 0,
	);
	my @snps = split(//,  $seq); #To iterate through the snps 

      foreach (@snps){ 
        if (my $tempkey = $codeHash{uc($_)}){
          foreach my $nuc(@{$tempkey}){ 
            $temp{$nuc}++; 
          } 
        } 
      }
	#Get frequency of ambiguities
	$temp{"N"} = ($seq =~ tr/NVHDB/NVHDB/ )*2;
	
	return(\%temp);
}

#Subroutine to get number of each nucleotide in an alignment
#Assumed HAPLOID (e.g. all het iupac codes are treated as ambiguity)
sub parseHapSeq{
	
	my $seq = shift(@_);
	#Assumes diploidy. iupac codes for 3 nuc ambiguity are treated as Ns (see sub percentn)
	my %codeHash =( 
		"A" => ["A"], 
		"T" => ["T"], 
		"G" => ["G"], 
		"C" => ["C"], 
		"-" => ["-"],
	); 
	#hash to store snp frequencies
	my %temp =(
		"A" => 0,
		"G"	=> 0,
		"C"	=> 0,
		"T"	=> 0,
		"N"	=> 0,
		"-"	=> 0,
	);
	my @snps = split(//,  $seq); #To iterate through the snps 

      foreach (@snps){ 
        if (my $tempkey = $codeHash{uc($_)}){
          foreach my $nuc(@{$tempkey}){ 
            $temp{$nuc}++; 
          } 
        } 
      }
	#Get frequency of ambiguities
	$temp{"N"} = ($seq =~ tr/RYKMSWNVHDB/RYKMSWNVHDB/);
	
	return(\%temp);
}

# subroutine to remove columns from an alignment, given the alignment contained in a hash and an array of positions in each value to be removed

sub removeColumns{

  my( $hashref, $remove ) = @_;

  my @blacklist = uniq($remove);
  
  # replace columns to be removed with a special character
  foreach my $key( sort keys %{ $hashref } ){
    foreach my $item( @blacklist ){
      substr(${$hashref}{$key}, $item, 1) = 'z';
    }
  }
  
  # replace the special characters with nothing
  foreach my $key( sort keys %{ $hashref } ){
    ${$hashref}{$key} =~ s/z//g;
  }
}


sub uniq {
	my @arr = @{$_[0]};
    my %seen;
    grep !$seen{$_}++, @arr;
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


# subroutine to report content of missing data in matrix

sub reportMatrixContent{
	
	my $hashref = shift;
	my $total = 0;
	my $missing = 0;
	my $gaps = 0;
	my $ambigs = 0;
	my $hets = 0;
	my $len;
	my $ind = scalar keys %{$hashref};
	
	foreach my $key (keys %{$hashref}){
	  $len = length(${$hashref}{$key});
      $total += length(${$hashref}{$key});   
	  $missing += (${$hashref}{$key} =~ tr/NVHDB/NVHDB/ );
	  $ambigs += (${$hashref}{$key} =~ tr/RYKMSWVHDB/RYKMSWVHDB/ );
	  $hets += (${$hashref}{$key} =~ tr/RYKMSW/RYKMSW/ );
	  $gaps += (${$hashref}{$key} =~ tr/-/-/ );
	}
	
	$total <= 0 and die "\nTotal nucleotides in data matrix less than or equal to zero. Something went wrong.\n\n";
	my $pMiss = $missing / $total * 100;
	my $pGap = $gaps / $total * 100;
	my $pAmbig = $ambigs / $total * 100;
	my $pHet = $hets / $total * 100;
	
	print "----------------REPORT----------------\n";
	print "Total remaining individuals = " . $ind . "\n";
	print "Total remaining nucleotide columns = " . $len . "\n";
	
	print "Total percent missing data in matrix = " ;
	printf("%.2f", $pMiss);
	print "%\nTotal percent ambiguities (excluding Ns) = " ;
	printf("%.2f", $pAmbig);
	print "%\nTotal percent heterozygous sites  = " ;
	printf("%.2f", $pHet);
	print "%\nTotal percent gaps in matrix = ";
	printf("%.2f", $pGap);
	print "%\n\n";
	
}










