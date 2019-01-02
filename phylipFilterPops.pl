#! /usr/bin/perl

# Contributions by Tyler K. Chafin, Steven M. Mussmann, Max R. Bangs
# Contact: tkchafin@uark.edu

use strict;
use warnings;
use Getopt::Std;
use List::Util qw(shuffle);

#Die if no arguments given
if( scalar( @ARGV ) == 0 ){
  &help;
  die "No options given\n\n";
}

#Parse arguments
my %opts;
getopts( 'p:i:o:hn:N:gPxrOm:sr:w:R:P', \%opts );

# kill if help option is true
if( $opts{h} ){
  &help;
  die "Help menu\n\n";
}

#get options
my ($map, $phy, $out, $threshold, $globalThresh, $gapFalse, $minPop, $skipMiss, $randSample, $weightSample, $alleleSample, $poplabels) = &parseArgs(\%opts);


#sanity check
if ($weightSample && $randSample || $weightSample && $alleleSample || $alleleSample && $randSample){
	die("-r, -R, and -w options cannot be combined\n")
}

if ($randSample){
	print("Randomly sampling $randSample samples per population\n");
}
if ($weightSample){
	print("Sampling $weightSample samples of highest coverage per population\n");
}
if ($alleleSample){
	print("Randomly sampling $alleleSample alleles per population (not working yet)\n");
}

# hash of loci with too much missing data
my %blacklist;

#parse popmap file into a hash with ind as key and popID as value
my $assignRef = &parsePopmap($map);

#parse phylip file into a hash with ind as key and array of seqs as value
my ($allRef, $ntax, $nchar) = &parsePhylip($phy);

#Print argument report
print "\nPopmap file is: $map\n";
print "Phylip file is: $phy\n";
print "Total taxa in phylip file: $ntax\n";
print "Total characters in phylip data matrix: $nchar\n\n";

#Get pop alignments only with ind as key and array of seqs as value
my @vals = values %$assignRef;
my @pops = uniq(@vals);

my $popAligns = &getSepMultPops($assignRef, $allRef, \@pops);

#remove pops with too little data
foreach my $p (keys %{$popAligns}){
	my $count=0;
	foreach my $ind (keys %{$popAligns->{$p}}){
		$count++;
	}
	if ($count < $minPop){
		print "Deleting population $p: Less than $minPop samples (N=$count).\n";
		delete $popAligns->{$p};
	}
}

# foreach my $p (keys %{$popAligns}){
# 	print $p, "\n";
# }
if ($skipMiss == 0){
	if ($gapFalse == 1){
		&getBlacklistSep($popAligns, $threshold, $globalThresh, \%blacklist);
	}else{
		&getBlacklistGap($popAligns, $threshold, $globalThresh, \%blacklist);
	}
	my $nFail = (keys %blacklist);

	print($nFail ," loci had greater than ",$threshold, " missing data. Removing them.\n");
}else{
	print "Skipping calculation of missing data.\n";
}

print("Writing new PHYLIP file <$out>\n");
open (PHY, "> $out");
my $locnum = 0;
my $indnum = 0;
for (my $loc = 0; $loc < $nchar; $loc++){
	if(!exists $blacklist{$loc+1}){
		$locnum++;
	}
}

foreach my $pop (keys %{$popAligns}){
	foreach my $ind (keys %{$popAligns->{$pop}}){
		$indnum++;
	}
}

print PHY $indnum, " ", $locnum, "\n";

#Get missing data per sample if needed
my %sampleCoverage;
if ($weightSample){
	%sampleCoverage = &getCoverage($popAligns);
	# foreach my $ind (keys %sampleCoverage){
	# 	print $ind, ": ", $sampleCoverage{$ind}, "\n";
	# }
	# exit()
}


#print data
foreach my $pop (keys %{$popAligns}){
	my $popCounter = 1;

	#if random sample, randomly order keys to continue, and only choose X number
	#if allele sample, grab all alleles and sample up to X non-N alleles
	#if weight sample, order keys by completeness of sampling (%sampCoverage)

	my @keys;
	my $max;
	if ($weightSample){
		@keys = sort { $sampleCoverage{$a} <=> $sampleCoverage{$b} } keys %{$popAligns->{$pop}};
		$max = $weightSample;
	}elsif ($randSample){
		@keys = shuffle (keys %{$popAligns->{$pop}});
		$max = $randSample;
	}else{
		@keys = keys %{$popAligns->{$pop}};
		$max = scalar(@keys);
	}

	#If not rand allele sampling, then sample based on the above order
	if (!$alleleSample){
		my $sampled = 0;
		foreach my $ind (@keys){
			$sampled++;
			#print $ind, "\n";
			my $name;
			if ($poplabels){
				if ($weightSample || $randSample){
					if ($max > 1){
						$name = $pop + "_" + $sampled;
					}else{
						$name = $pop;
					}
				}
			}else{
				$name = $ind;
			}
			print PHY $name, "\t";
			for (my $l = 0; $l < $nchar; $l++){
				if(!exists $blacklist{$l+1}){
					#check if samples should be relabeled
					#if -r, -R, or -w equal 1, don't append number
					print PHY ${$popAligns->{$pop}->{$ind}}[$l];
				}
			}
			print PHY "\n";
			if ($sampled >= $max){
				last;
			}
		}
		$popCounter++; #increment count per pop
	}
}

close PHY;




exit 0;

 ########################### SUBROUTINES ###############################

 sub help{

	print "\nphylipFilterPops.pl by Tyler Chafin, w/ contributions by Steve Mussmann and Max Bangs\n";
	print "\nThis script filters rows and columns from a phylip file based on population assignments\n";
	print "A population map should be given in a tab-delimited file, formatted as:\n";
	print "\n\tSampleName\tPopID\n\n";
	print "Where PopID can be a string or integer, and SampleName must exactly match a corresponsing line in the phylip file. Any samples not in the popmap will not be included in the output files.\n\n";
	print "Options:\n";
	print "\t-p	: Path to popmap file (tab-delimited)\n";
	print "\t-i	: Path to input file (phylip)\n";
	print "\t-m	: Minimum number of samples to keep a population\n";
	print "\t-n	: Proportion missing data allowed per population per SNP (default=0.5)\n";
	print "\t-N	: Proportion of globally missing data allowed per SNP (default=0.5)\n";
	print "\t-g	: Toggle on to TURN OFF default behavior of treating gaps as missing data\n";
	print "\t-s	: Skip calculating missing data, and just drop populations with too few inds\n";
	print "\t-r	: Randomly sample up to <x> individuals from each population\n";
	print "\t-w	: Sample up to <x> individuals from each population with the least missing data\n";
	print "\t-R	: Randomly sample up to <x> ALLELES from each population (per column, for SNPs)\n";
	print "\t-P	: Print samples using pop names (numbered as \"_x\" if more than 1 sample)\n";
	print "\t-o	: Output file name\n";
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
	my $threshold  = $opts{n} || 0.5;
	my $gapFalse = 0;
	my $phyNew = 0;
	my $skip = 0;
	$opts{s} and $skip = 1;
	my $minPop = $opts{m} || 0;
	$opts{g} and $gapFalse = 1;
	my $globalThresh = $opts{N} || 0.5;
  my $out = $opts{o} || "out.phy";
	my $randSample = $opts{r} || 0;
	my $weightSample = $opts{w} || 0;
	my $alleleSample = $opts{R} || 0;
	my $poplabels = 0;
	$opts{P} and $poplabels = 1;
  #return
  return ($map, $phy, $out, $threshold, $globalThresh, $gapFalse, $minPop, $skip, $randSample, $weightSample, $alleleSample, $poplabels);
}

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
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
sub getSepMultPops{

	my $popRef 		= $_[0];
	my $seqRef 		= $_[1];
	my $toGetRef 	= $_[2];

	my %pop; #hash of hashes

	foreach my $id (@{$toGetRef}){
		#print $id, "\n";
		my %local;
		$pop{$id} = \%local;
	}

	foreach my $id (@{$toGetRef}){
		foreach my $key (keys %{$popRef}){
			#If pop ID matches, get sequence
			if ($popRef->{$key} eq $id){
				if (exists $seqRef->{$key}){
					$pop{$id}{$key} = $seqRef->{$key};
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



sub calcMissingSep{

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

sub getCoverage{

	my $popsRef = $_[0];
	my %coverage;

	foreach my $pop (keys %{$popsRef}){
		foreach my $ind (keys %{$popsRef->{$pop}}){
			my $totalcount = 0;
			my $ncount = 0;
			foreach my $locus ( @{$popsRef->{$pop}->{$ind}} ){
				if ($locus eq "N"){
					$ncount++;
				}
				$totalcount++;
			}
			$coverage{$ind} = ($ncount/$totalcount);
		}
	}
	return(%coverage);
}

sub getBlacklistSep{

	my( $popsRef, $thresh, $globalThresh, $blacklistref ) = @_;

	my $globalInds = 0;
	my %globalCount;

	#Check loci in each population
	foreach my $pop (keys %{$popsRef}){
		my %ncount;
		my $inds;
		foreach my $ind (keys %{$popsRef->{$pop}}){
			$inds++;
			$globalInds++;
			my $counter = 0;
			foreach my $locus ( @{$popsRef->{$pop}->{$ind}} ){
				$counter++;
				$ncount{$counter} = 0 unless exists $ncount{$counter};
				$globalCount{$counter} = 0 unless exists $globalCount{$counter};
				if ($locus eq "N"){
					$ncount{$counter}++;
					$globalCount{$counter}++;
				}
			}
		}
		#blacklist any loci with too high n proportion (in THIS pop)
		foreach my $loc(sort keys %ncount){
			if (($ncount{$loc} / $inds) > $threshold){
				$$blacklistref{$loc} = ($ncount{$loc} / $inds) unless exists $$blacklistref{$loc};
				#print("Locus:",$loc," - Failed (pop", $pop,"): ",($ncount{$loc} / $inds), "\n");
			}
			#print("Nprop:",($ncount{$loc} / $inds), "\n");
		}
	}

	#Blacklist loci with globally too high missing data
	foreach my $loc(sort keys %globalCount){
		if (($globalCount{$loc} / $globalInds) > $threshold){
			$$blacklistref{$loc} = ($globalCount{$loc} / $globalInds) unless exists $$blacklistref{$loc};
			#print("Locus:",$loc," - Failed (global): ",($globalCount{$loc} / $globalInds), "\n");
		}
	}
}

sub getBlacklistGap{

	my( $popsRef, $thresh, $globalThresh, $blacklistref ) = @_;

	my $globalInds = 0;
	my %globalCount;

	#Check loci in each population
	foreach my $pop (keys %{$popsRef}){
		my %ncount;
		my $inds;
		foreach my $ind (keys %{$popsRef->{$pop}}){
			#print("Ind: ", $ind, "\n");
			$inds++;
			$globalInds++;
			my $counter = 0;
			foreach my $locus ( @{$popsRef->{$pop}->{$ind}} ){
				$counter++;
				$ncount{$counter} = 0 unless exists $ncount{$counter};
				$globalCount{$counter} = 0 unless exists $globalCount{$counter};
				if ($locus eq "N" || $locus eq "-"){
					$ncount{$counter}++;
					$globalCount{$counter}++;
				}
			}
		}
		#blacklist any loci with too high n proportion (in THIS pop)
		foreach my $loc(sort keys %ncount){
			if (($ncount{$loc} / $inds) > $threshold){
				$$blacklistref{$loc} = ($ncount{$loc} / $inds) unless exists $$blacklistref{$loc};
				#print("Locus:",$loc," - Failed (pop", $pop,"): ",($ncount{$loc} / $inds), "\n");
			}
			#print("Nprop:",($ncount{$loc} / $inds), "\n");
		}
	}

	#Blacklist loci with globally too high missing data
	foreach my $loc(sort keys %globalCount){
		if (($globalCount{$loc} / $globalInds) > $threshold){
			$$blacklistref{$loc} = ($globalCount{$loc} / $globalInds) unless exists $$blacklistref{$loc};
			#print("Locus:",$loc," - Failed (global): ",($globalCount{$loc} / $globalInds), "\n");
		}
	}

}
