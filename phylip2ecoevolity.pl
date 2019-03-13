#! /usr/bin/perl

# Contributions by Tyler K. Chafin, Steven M. Mussmann, Max R. Bangs
# Contact: tkchafin@uark.edu

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

#Die if no arguments given
if( scalar( @ARGV ) == 0 ){
  &help;
  die "No options given\n\n";
}

#Parse arguments
my %opts;
getopts( 'p:i:1:2:sn:N:gPo:h', \%opts );

# kill if help option is true
if( $opts{h} ){
  &help;
  die "Help menu\n\n";
}

#get options
my ($map, $phy, $pop1, $pop2, $threshold, $globalThresh, $gapFalse, $phyNew, $out, $same) = &parseArgs(\%opts);

#Extract pops into an array
my @pop1list = split(/\+/,$pop1);
my @pop2list = split(/\+/,$pop2);

my $a1 = join(', ', @pop1list);
my $a2 = join(', ', @pop2list);

# hash for loci with too much missing data
my %blacklist;

#parse popmap file into a hash with ind as key and popID as value
my $assignRef = &parsePopmap($map);

#parse phylip file into a hash with ind as key and array of seqs as value
my ($allRef, $ntax, $nchar) = &parsePhylip($phy);

#Print argument report
print "\nPopmap file is: $map\n";
print "Phylip file is: $phy\n";
if (scalar @pop2list > 0){
	print "Selected populations (1 pop model): $a1\n";
}else{
	print "Selected populations (2 pop model): $a1 and $a2\n";
}
print "Total taxa in phylip file: $ntax\n";
print "Total characters in phylip data matrix: $nchar\n\n";

#Get pop alignments only with ind as key and array of seqs as value
my $pop1Ref = &getMultPops($assignRef, $allRef, \@pop1list);
my $pop2Ref = "";
if (scalar @pop2list > 0){
	$pop2Ref = &getMultPops($assignRef, $allRef, \@pop2list);
}
# calculate missing data in admixed population
#&calcMissing($popaRef, \%blacklist);
if ($gapFalse == 1){
	&getBlacklist($threshold, scalar(@pop2list), \%blacklist, $pop1Ref, $pop2Ref);
}else{
	&getBlacklistGap($threshold, scalar(@pop2list), \%blacklist, $pop1Ref, $pop2Ref);
}
my $nFail = (keys %blacklist);

print($nFail ," loci had greater than ",$threshold, " missing data. Removing them.\n");
#print Dumper(\%blacklist);

#Check if pops contain data
my $num1 = keys %{$pop1Ref};
my $num2 = keys %{$pop2Ref};

if ($num1< 1){
  die "Error: No individuals were found <$a1>!\n\n";
}else{
  print "Found <$num1> individuals in selected pops\n";
}
if ($num2< 1){
  die "Error: No individuals were found <$a2>!\n\n";
}else{
  print "Found <$num2> individuals in selected pops\n";
}

my $indnum = 0;
my $locnum = 0;
for (my $loc = 0; $loc < $nchar; $loc++){
	if(!exists $blacklist{$loc+1}){
		$locnum++;
	}
}
foreach my $ind (keys %{$pop1Ref}){
	$indnum++;
}
foreach my $ind (keys %{$pop2Ref}){
	$indnum++;
}

if ($phyNew){
	print("Writing new PHYLIP file <out.phy>\n");
	open (PHY, "> out.phy");
	print PHY $indnum, " ", $locnum, "\n";

	#print data for P1
	foreach my $ind (sort keys %{$pop1Ref}){
		#print $ind, "\n";
		print PHY "Pop1_",$ind, "\t";
		for (my $l = 0; $l < $nchar; $l++){
			if(!exists $blacklist{$l+1}){
				print PHY ${$pop1Ref->{$ind}}->[$l];
			}
		}
		print PHY "\n";
	}
	if (scalar @pop2list > 0){
		#print data for P1
		foreach my $ind (sort keys %{$pop1Ref}){
			#print $ind, "\n";
			print PHY "Pop1_",$ind, "\t";
			for (my $l = 0; $l < $nchar; $l++){
				if(!exists $blacklist{$l+1}){
					print PHY ${$pop1Ref->{$ind}}->[$l];
				}
			}
			print PHY "\n";
		}
	}
	close PHY;
}


exit 0;

 ########################### SUBROUTINES ###############################

 sub help{

	print "\nphylip2ecoevolity.pl\n";
	print "\nThis script converts from phylip format to the NEXUS input for Ecoevolity (Oaks, 2019)\n";
	print "A population map should be given in a tab-delimited file, formatted as:\n";
	print "\n\tSampleName\tPopID\n\n";
	print "Where PopID can be a string or integer, and SampleName must exactly match a corresponsing line in the phylip file. Any samples not in the popmap will not be included in the output files.\n\n";
	print "Options:\n";
	print "\t-p	: Path to popmap file (tab-delimited)\n";
	print "\t-1	: Identifier for populations to include in output (include multiple as: pop1+pop2)\n";
	print "\t-2	: Identifier for populations to include in output, if using a 2-pop comparison\n";
	print "\t-s	: Toggle on to use the same population identifier for all output individuals.\n";
	print "\t-i	: Path to input file (phylip)\n";
	print "\t-n	: Proportion missing data allowed per pop per column (default=0.1)\n";
	print "\t-n	: Proportion missing data allowed globally per column (default=0.1)\n";
	print "\t-g	: Toggle on to TURN OFF default behavior of treating gaps as missing data\n";
	print "\t-P	: Toggle on to output a new phylip file with the filtered data. [default=off]\n";
	print "\t-o	: Output prefix [default = \"ee_out\"]\n";
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
  my $pops1  = $opts{1} or die "\nNo population specified.\n\n";
	my $pops2 = "";
	$pops2  = $opts{2};
	my $threshold  = $opts{n} || 0.1;
	my $globalThresh = $opts{N} || 0.1;
	my $gapFalse = 0;
	my $phyNew = 0;
	$opts{g} and $gapFalse = 1;
	$opts{P} and $phyNew = 1;
  my $out = $opts{o} || "ee_out";
	my $same = 0;
	$opts{s} and $same = 1;
  #return
  return ($map, $phy, $pops1, $pops2, $threshold, $globalThresh, $gapFalse, $phyNew, $out, $same);
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
					print "Warning: Sample <$key> was not found in sequence file. Skipping.\n";
				}
			}
		}
	}
	return(\%pop);
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
	my @blacklist;

	#To track fixed alleles in each pop
	my $alleles1 = parseColumn($p1, $thresholdN, $thresholdG, \@blacklist);
	my $alleles2 = parseColumn($p2, $thresholdN, $thresholdG, \@blacklist);

	#Make sure both pops have same number of columns
	if ((scalar(@{$alleles1})) != (scalar(@{$alleles1}))){
		die "\nError: Y ur populations have not same sequence leNGTH???\n\n";
	}else{
		#Only keep loci which are differentially fixed
		#Make sure to check anything fixed in pop1 is different
		#from fixed in pop2
		for(my $i=0; $i < scalar(@{$alleles1}); $i++){
			my $check1 = $alleles1->[$i] =~ tr/NV-/NV-/;
			my $check2 = $alleles2->[$i] =~ tr/NV-/NV-/;
			#If either pop was variable, or fixed for gaps or Ns
			if ($check1 > 0 || $check2 > 0){
				next;
			}else{
				#If both fixed for same allele
				if ($alleles1->[$i] eq $alleles2->[$i]){
					push(@blacklist, $i);
					next;
				}
			}
		}
	}
	return(\@blacklist);
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

sub calcMissing{

	my( $hashref, $blacklistref ) = @_;

	foreach my $ind( sort keys %$hashref ){
		my $counter = 0;
		foreach my $locus( @${$$hashref{$ind}} ){
			$counter++;
			if($locus eq "N"){
				$$blacklistref{$counter}++;
			}else{
				$$blacklistref{$counter}+=0;
			}
		}
	}
}

sub getBlacklistGap{

	my( $thresh, $p2, $blacklistref, $p1ref, $p2ref) = @_;

	my $globalInds = 0;
	my %globalCount;
	
	#check loci in admixed
	my %ncount;
	my $inds1;
	foreach my $ind( sort keys %$p1ref ){
		$inds1++;
		$globalInds++;
		my $counter = 0;
		foreach my $locus( @${$$p1ref{$ind}} ){
			$counter++;
			$ncount{$counter} = 0 unless exists $ncount{$counter};
			$globalCount{$counter} = 0 unless exists $globalCount{$counter};
			if($locus eq "N"){
				$ncount{$counter}++;
				$globalCount{$counter}++;
			}
		}
	}
	#blacklist any loci with too high n proportion
	foreach my $loc(sort keys %ncount){
		if (($ncount{$loc} / $inds1) > $threshold){
			$$blacklistref{$loc} = ($ncount{$loc} / $inds1) unless exists $$blacklistref{$loc};
			#print("Locus:",$loc," - Missing data: ",($nAcount{$loc} / $admixinds), "\n");
		}
	}
	
	if ($p2 > 0){
		my %ncount2;
		my $inds2;
		foreach my $ind( sort keys %$p2ref ){
			$inds2++;
			$globalInds++;
			my $counter = 0;
			foreach my $locus( @${$$p2ref{$ind}} ){
				$counter++;
				$ncount2{$counter} = 0 unless exists $ncount2{$counter};
				$globalCount{$counter} = 0 unless exists $globalCount{$counter};
				if($locus eq "N"){
					$ncount2{$counter}++;
					$globalCount{$counter}++;
				}
			}
		}
		#blacklist any loci with too high n proportion
		foreach my $loc(sort keys %ncount2){
			if (($ncount2{$loc} / $inds2) > $threshold){
				$$blacklistref{$loc} = ($ncount2{$loc} / $inds2) unless exists $$blacklistref{$loc};
				#print("Locus:",$loc," - Missing data: ",($nAcount{$loc} / $admixinds), "\n");
			}
		}
		#Blacklist loci with globally too high missing data
		foreach my $loc(sort keys %globalCount){
			if (($globalCount{$loc} / $globalInds) > $threshold){
				$$blacklistref{$loc} = ($globalCount{$loc} / $globalInds) unless exists $$blacklistref{$loc};
				#print("Locus:",$loc," - Missing data: ",($nAcount{$loc} / $admixinds), "\n");
			}
		}
	}
}

sub getBlacklistGap{

	my( $thresh, $p2, $blacklistref, $p1ref, $p2ref) = @_;

	my $globalInds = 0;
	my %globalCount;
	
	#check loci in admixed
	my %ncount;
	my $inds1;
	foreach my $ind( sort keys %$p1ref ){
		$inds1++;
		$globalInds++;
		my $counter = 0;
		foreach my $locus( @${$$p1ref{$ind}} ){
			$counter++;
			$ncount{$counter} = 0 unless exists $ncount{$counter};
			$globalCount{$counter} = 0 unless exists $globalCount{$counter};
			if($locus eq "N" || $locus eq "-"){
				$ncount{$counter}++;
				$globalCount{$counter}++;
			}
		}
	}
	#blacklist any loci with too high n proportion
	foreach my $loc(sort keys %ncount){
		if (($ncount{$loc} / $inds1) > $threshold){
			$$blacklistref{$loc} = ($ncount{$loc} / $inds1) unless exists $$blacklistref{$loc};
			#print("Locus:",$loc," - Missing data: ",($nAcount{$loc} / $admixinds), "\n");
		}
	}
	
	if ($p2 > 0){
		my %ncount2;
		my $inds2;
		foreach my $ind( sort keys %$p2ref ){
			$inds2++;
			$globalInds++;
			my $counter = 0;
			foreach my $locus( @${$$p2ref{$ind}} ){
				$counter++;
				$ncount2{$counter} = 0 unless exists $ncount2{$counter};
				$globalCount{$counter} = 0 unless exists $globalCount{$counter};
				if($locus eq "N" || $locus eq "-"){
					$ncount2{$counter}++;
					$globalCount{$counter}++;
				}
			}
		}
		#blacklist any loci with too high n proportion
		foreach my $loc(sort keys %ncount2){
			if (($ncount2{$loc} / $inds2) > $threshold){
				$$blacklistref{$loc} = ($ncount2{$loc} / $inds2) unless exists $$blacklistref{$loc};
				#print("Locus:",$loc," - Missing data: ",($nAcount{$loc} / $admixinds), "\n");
			}
		}
		#Blacklist loci with globally too high missing data
		foreach my $loc(sort keys %globalCount){
			if (($globalCount{$loc} / $globalInds) > $threshold){
				$$blacklistref{$loc} = ($globalCount{$loc} / $globalInds) unless exists $$blacklistref{$loc};
				#print("Locus:",$loc," - Missing data: ",($nAcount{$loc} / $admixinds), "\n");
			}
		}
	}
}
