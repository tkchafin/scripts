#! /usr/bin/perl

# Contributions by Max Bangs, Steve Mussmann, and Tyler Chafin
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
getopts( 'p:i:1:2:a:o:hn:N:gPxr:', \%opts );

# kill if help option is true
if( $opts{h} ){
  &help;
  die "Help menu\n\n";
}

#get options
my ($map, $phy, $p1, $p2, $ad, $out, $threshold, $globalThresh, $gapFalse, $phyNew, $onlyPhy, $rand) = &parseArgs(\%opts);

#Extract pops into an array
my @pop1 = split(/\+/,$p1);
my @pop2 = split(/\+/,$p2);
my @popA = split(/\+/,$ad);

$p1 = join(', ',@pop1);
$p2 = join(', ',@pop2);
$ad = join(', ', @popA);

# hash of loci with too much missing data
my %blacklist;

#parse popmap file into a hash with ind as key and popID as value
my $assignRef = &parsePopmap($map);

#parse phylip file into a hash with ind as key and array of seqs as value
my ($allRef, $ntax, $nchar) = &parsePhylip($phy);

#Print argument report
print "\nPopmap file is: $map\n";
print "Phylip file is: $phy\n";
print "Population 1 is: $p1\n";
print "Population 2 is: $p2\n";
print "Admixed population is: $ad\n";
print "Total taxa in phylip file: $ntax\n";
print "Total characters in phylip data matrix: $nchar\n\n";

#Get pop alignments only with ind as key and array of seqs as value
my $pop1Ref = &getMultPops($assignRef, $allRef, \@pop1);
my $pop2Ref = &getMultPops($assignRef, $allRef, \@pop2);
my $popaRef = &getMultPops($assignRef, $allRef, \@popA);

# calculate missing data in admixed population
#&calcMissing($popaRef, \%blacklist);
if ($gapFalse == 1){
	&getBlacklist($threshold, $globalThresh, $pop1Ref, $pop2Ref, $popaRef, \%blacklist);
}else{
	&getBlacklistGap($threshold, $globalThresh, $pop1Ref, $pop2Ref, $popaRef, \%blacklist);
}
my $nFail = (keys %blacklist);

print($nFail ," loci had greater than ",$threshold, " missing data. Removing them.\n");


#Random sample from remaining loci
if ($rand){
	print "Randomly sampling ",$rand," SNPs...\n";
	&randomSample($nchar, \%blacklist, $rand); #Function adds unsampled loci to the blacklist
	my $rem = scalar keys %blacklist;
	print "After random sampling and filtering, a total of ", $rem, " SNPs will be removed!\n";
}

#print Dumper(\%blacklist);

#Check if pops contain data
my $num1 = keys %{$pop1Ref};
my $num2 = keys %{$pop2Ref};
my $num3 = keys %{$popaRef};

if ($num1< 1){
  die "Error: No individuals for pop ID <$p1> were found!\n\n";
}elsif ($num2 < 1){
  die "Error: No individuals for pop ID <$p2> were found!\n\n";
}elsif ($num3 < 1){
  die "Error: No individuals for pop ID <$ad> were found!\n\n";
}else{
  print "Found <$num1> individuals in population 1\n";
  print "Found <$num2> individuals in population 2\n";
  print "Found <$num3> individuals in admixed population\n\n";
}

#Open filstreams
my $locnum = 0;
my $indnum = 0;
if ($onlyPhy != 1){
	open(NEWHYB, "> newhyb_dat.txt");

	#Write header
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
	foreach my $ind (keys %{$popaRef}){
		$indnum++;
	}
	print NEWHYB "NumIndivs ",$indnum, "\n";
	print NEWHYB "NumLoci ",$locnum, "\n";
	print NEWHYB "Digits 1\nFormat Lumped\n";
	print NEWHYB "LocusNames\t";
	for my $nloci (1 .. $nchar){
		if(!exists $blacklist{$nloci}){
			print NEWHYB "l", "$nloci\t";
		}
	}
	print NEWHYB "\n";

	#Print P1 pure individuals first
	open(IND, "> ind_order.txt");
	open(POPS, ">pop_order.txt");

	print"Outputting individuals to NewHybrids file (in THIS order):\n";
	my $count = 1;
	foreach my $ind (sort keys %{$pop1Ref}){
		print $ind, "\n";
		print IND $ind, "\n";
		print POPS $assignRef->{$ind}, "\n";
		print NEWHYB $count, "\tz0\t";
		for (my $loc = 0; $loc < $nchar; $loc++){
			if(!exists $blacklist{$loc+1}){
				my $nuc = uc ${$pop1Ref->{$ind}}->[$loc];
				my $new;
				if ($nuc eq "A") {$new = "11"}
				elsif ($nuc eq "T") {$new = "22"}
				elsif ($nuc eq "G") {$new = "33"}
				elsif ($nuc eq "C") {$new = "44"}
				elsif ($nuc eq "W") {$new = "12"}
				elsif ($nuc eq "R") {$new = "13"}
				elsif ($nuc eq "M") {$new = "14"}
				elsif ($nuc eq "K") {$new = "23"}
				elsif ($nuc eq "Y") {$new = "24"}
				elsif ($nuc eq "S") {$new = "34"}
				elsif ($nuc eq "N") {$new = "0"}
				else {$new = "0"}
				print NEWHYB $new, "\t";
			}
		}
		print NEWHYB "\n";
		$count++;
	}
	#Print P2 pures
	foreach my $ind (sort keys %{$pop2Ref}){
		print $ind, "\n";
		print IND $ind, "\n";
		print POPS $assignRef->{$ind}, "\n";
		print NEWHYB $count, "\tz1\t";
		for (my $loc = 0; $loc < $nchar; $loc++){
			if(!exists $blacklist{$loc+1}){
				my $nuc = uc ${$pop2Ref->{$ind}}->[$loc];
				my $new;
				if ($nuc eq "A") {$new = "11"}
				elsif ($nuc eq "T") {$new = "22"}
				elsif ($nuc eq "G") {$new = "33"}
				elsif ($nuc eq "C") {$new = "44"}
				elsif ($nuc eq "W") {$new = "12"}
				elsif ($nuc eq "R") {$new = "13"}
				elsif ($nuc eq "M") {$new = "14"}
				elsif ($nuc eq "K") {$new = "23"}
				elsif ($nuc eq "Y") {$new = "24"}
				elsif ($nuc eq "S") {$new = "34"}
				elsif ($nuc eq "N") {$new = "0"}
				else {$new = "0"}
				print NEWHYB $new, "\t";
			}
		}
		print NEWHYB "\n";
		$count++;
	}
	#Print admixed individuals
	foreach my $ind (sort keys %{$popaRef}){
		print $ind, "\n";
		print IND $ind, "\n";
		print POPS $assignRef->{$ind}, "\n";
		print NEWHYB $count, "\t\t";
		for (my $loc = 0; $loc < $nchar; $loc++){
			if(!exists $blacklist{$loc+1}){
				my $nuc = uc ${$popaRef->{$ind}}->[$loc];
				my $new;
				if ($nuc eq "A") {$new = "11"}
				elsif ($nuc eq "T") {$new = "22"}
				elsif ($nuc eq "G") {$new = "33"}
				elsif ($nuc eq "C") {$new = "44"}
				elsif ($nuc eq "W") {$new = "12"}
				elsif ($nuc eq "R") {$new = "13"}
				elsif ($nuc eq "M") {$new = "14"}
				elsif ($nuc eq "K") {$new = "23"}
				elsif ($nuc eq "Y") {$new = "24"}
				elsif ($nuc eq "S") {$new = "34"}
				elsif ($nuc eq "N") {$new = "0"}
				else {$new = "0"}
				print NEWHYB $new, "\t";
			}
		}
		print NEWHYB "\n";
		$count++;
	}
	print "Done writing NEWHYB file newhyb_dat.txt>\n\n";
	close NEWHYB;
	close IND;
	close POP;
	print "Wrote <ind_order.txt> for order of individuals.\n";
	print "Wrore <pop_order.txt> reporting Pop IDs of individuals. Use this with newhybs2distruct.py\n";
}

if ($phyNew or $onlyPhy){
	print("Writing new PHYLIP file <out.phy>\n");
	open (PHY, "> out.phy");
	my $locnum = 0;
	my $indnum = 0;
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
	foreach my $ind (keys %{$popaRef}){
		$indnum++;
	}
	print PHY $indnum, " ", $locnum, "\n";

	#print data for P1
	foreach my $ind (sort keys %{$pop1Ref}){
		#print $ind, "\n";
		print PHY $ind, "\t";
		for (my $l = 0; $l < $nchar; $l++){
			if(!exists $blacklist{$l+1}){
				print PHY ${$pop1Ref->{$ind}}->[$l];
			}
		}
		print PHY "\n";
	}
	#exit;
	foreach my $ind (sort keys %{$pop2Ref}){
		print PHY $ind, "\t";
		for (my $loc = 0; $loc < $nchar; $loc++){
			#print $loc, " ";
			if(!exists $blacklist{$loc+1}){
				print PHY ${$pop2Ref->{$ind}}->[$loc];
			}
		}
		print PHY "\n";
	}
	foreach my $ind (sort keys %{$popaRef}){
		print PHY $ind, "\t";
		for (my $loc = 0; $loc < $nchar; $loc++){
			if(!exists $blacklist{$loc+1}){
				print PHY ${$popaRef->{$ind}}->[$loc];
			}
		}
		print PHY "\n";
	}
	close PHY;
}

exit 0;

 ########################### SUBROUTINES ###############################

 sub help{

	print "\nphylip2newhybrids.pl by Max Bangs, Steve Mussmann, and Tyler Chafin\n";
	print "\nThis script converts from phylip format to the input format for NewHybrids\n";
	print "A population map should be given in a tab-delimited file, formatted as:\n";
	print "\n\tSampleName\tPopID\n\n";
	print "Where PopID can be a string or integer, and SampleName must exactly match a corresponsing line in the phylip file. Any samples not in the popmap will not be included in the output files.\n\n";
	print "Options:\n";
	print "\t-p	: Path to popmap file (tab-delimited)\n";
	print "\t-1	: Identifier for population 1 (include multiple as: pop1+pop2)\n";
	print "\t-2	: Identifier for population 2 (include multiple as: pop1+pop2)\n";
	print "\t-a	: Identifier for admixed population(s) (include multiple as: pop1+pop2)\n";
	print "\t-i	: Path to input alignment file (phylip)\n";
	print "\t-n	: Proportion missing data allowed per population per SNP (default=0.5)\n";
	print "\t-N	: Proportion of globally missing data allowed per SNP (default=0.5)\n";
	print "\t-g	: Toggle on to TURN OFF default behavior of treating gaps as missing data\n";
	print "\t-P	: Toggle on to output a new phylip file with the filtered data. [default=off]\n";
	print "\t-r	: Randomly sample x number of SNPs post-filtering [default=off]";
	print "\t-x	: Toggle on to ONLY print a new phylip file (e.g. no NewHybrids file)";
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
  my $ad  = $opts{a} or die "\nNo admixed population specified.\n\n";
	my $threshold  = $opts{n} || 0.5;
	my $gapFalse = 0;
	my $phyNew = 0;
	$opts{g} and $gapFalse = 1;
	$opts{P} and $phyNew = 1;
	my $globalThresh = $opts{N} || 0.5;
	my $onlyPhy = 0;
	my $rand = $opts{r} || 0;
	$opts{x} and $onlyPhy = 1;
  my $out = $opts{o} || "out.phy";
  #return
  return ($map, $phy, $p1, $p2, $ad, $out, $threshold, $globalThresh, $gapFalse, $phyNew, $onlyPhy, $rand);
}

#################################################################################
#subroutine to randomly sample SNPs to keep in alignment
sub randomSample{

 my( $nchar, $blacklistref, $num ) = @_;

 my $total = $nchar - (scalar keys %$blacklistref);

 if ($num >= $total){
 print "\n[-r] invoked, but number to sample is greater than total SNPs remaining after filtering!\n\n";
	 return(1);
 }

 my @indices;
 for my $loc (1 .. $nchar){
	 if(!exists $blacklist{$loc}){
		 push(@indices, $loc);
	 }
 }

 my $tokill = (scalar @indices) - $num;

 shuffle(\@indices); #Shuffle indices


 my @remove = @indices[0..$tokill-1]; #Choose n of them for removal

 #print "Need to keep:",$num, "\n";
 #print "Total elements:",$total, "\n";
 #print "To delete:",$tokill, "\n";

 #  foreach my $k(@remove){
	#  print $k, ",";
 #  }
 #  print "\n";
 #
	# my %hash;
	# @hash{@remove}=();
	# print "Keeping: ";
 # foreach my $ind(@indices){
	#  if (!exists $hash{$ind}){
	# 	 print $ind, ",";
	#  }
 # }
 # print "\n";
 #my $rem = scalar @remove;
 #print("Removing:", $rem, "\n");

 foreach my $i (@remove){
	$$blacklistref{$i} = "randomfail";
 }
 return (0);

}

sub shuffle{
 my $array = shift;
 my $i = @$array;
 while (--$i){
	 my $j = int rand($i+1);
	 @$array[$i,$j] = @$array[$j,$i];
 }
}

#################################################################################



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
  open( OUT, '>', $out ) or die "Can't open $out!: $!\n\n";

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

sub getBlacklist{

	my( $thresh, $globalThresh, $p1ref, $p2ref, $admixref, $blacklistref ) = @_;

	my $globalInds = 0;
	my %globalCount;

	#Check loci in pop1
	my %ncount;
	my $p1inds;
	foreach my $ind( sort keys %$p1ref ){
		$p1inds++;
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
		if (($ncount{$loc} / $p1inds) > $threshold){
			$$blacklistref{$loc} = ($ncount{$loc} / $p1inds) unless exists $$blacklistref{$loc};
			#print("Locus:",$loc," - Missing data: ",($ncount{$loc} / $p1inds), "\n");
		}
		#print("Nprop:",($ncount{$loc} / $p1inds), "\n");
	}
	#check loci in pop2
	my %n2count;
	my $p2inds;
	foreach my $ind( sort keys %$p2ref ){
		$p2inds++;
		$globalInds++;
		my $counter = 0;
		foreach my $locus( @${$$p2ref{$ind}} ){
			$counter++;
			$n2count{$counter} = 0 unless exists $n2count{$counter};
			if($locus eq "N"){
				$n2count{$counter}++;
				$globalCount{$counter}++;
			}
		}
	}
	#blacklist any loci with too high n proportion
	foreach my $loc(sort keys %ncount){
		if (($n2count{$loc} / $p2inds) > $threshold){
			$$blacklistref{$loc} = ($n2count{$loc} / $p2inds) unless exists $$blacklistref{$loc};
			#print("Locus:",$loc," - Missing data: ",($n2count{$loc} / $p2inds), "\n");
		}
		#print("Nprop:",($n2count{$loc} / $p2inds), "\n");
	}
	#check loci in admixed
	my %nAcount;
	my $admixinds;
	foreach my $ind( sort keys %$admixref ){
		$admixinds++;
		$globalInds++;
		my $counter = 0;
		foreach my $locus( @${$$admixref{$ind}} ){
			$counter++;
			$nAcount{$counter} = 0 unless exists $nAcount{$counter};
			if($locus eq "N"){
				$nAcount{$counter}++;
				$globalCount{$counter}++;
			}
		}
	}
	#blacklist any loci with too high n proportion
	foreach my $loc(sort keys %ncount){
		if (($nAcount{$loc} / $admixinds) > $threshold){
			$$blacklistref{$loc} = ($nAcount{$loc} / $admixinds) unless exists $$blacklistref{$loc};
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

sub getBlacklistGap{

	my( $thresh, $globalThresh, $p1ref, $p2ref, $admixref, $blacklistref ) = @_;

	my $globalInds = 0;
	my %globalCount;

	#Check loci in pop1
	my %ncount;
	my $p1inds;
	foreach my $ind( sort keys %$p1ref ){
		$p1inds++;
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
		if (($ncount{$loc} / $p1inds) > $threshold){
			$$blacklistref{$loc} = ($ncount{$loc} / $p1inds) unless exists $$blacklistref{$loc};
			#print("Locus:",$loc," - Missing data: ",($ncount{$loc} / $p1inds), "\n");
		}
		#print("Nprop:",($ncount{$loc} / $p1inds), "\n");
	}
	#check loci in pop2
	my %n2count;
	my $p2inds;
	foreach my $ind( sort keys %$p2ref ){
		$p2inds++;
		$globalInds++;
		my $counter = 0;
		foreach my $locus( @${$$p2ref{$ind}} ){
			$counter++;
			$n2count{$counter} = 0 unless exists $n2count{$counter};
			if($locus eq "N" || $locus eq "-"){
				$n2count{$counter}++;
				$globalCount{$counter}++;
			}
		}
	}
	#blacklist any loci with too high n proportion
	foreach my $loc(sort keys %ncount){
		if (($n2count{$loc} / $p2inds) > $threshold){
			$$blacklistref{$loc} = ($n2count{$loc} / $p2inds) unless exists $$blacklistref{$loc};
			#print("Locus:",$loc," - Missing data: ",($n2count{$loc} / $p2inds), "\n");
		}
		#print("Nprop:",($n2count{$loc} / $p2inds), "\n");
	}
	#check loci in admixed
	my %nAcount;
	my $admixinds;
	foreach my $ind( sort keys %$admixref ){
		$admixinds++;
		$globalInds++;
		my $counter = 0;
		foreach my $locus( @${$$admixref{$ind}} ){
			$counter++;
			$nAcount{$counter} = 0 unless exists $nAcount{$counter};
			if($locus eq "N" || $locus eq "-"){
				$nAcount{$counter}++;
				$globalCount{$counter}++;
			}
		}
	}
	#blacklist any loci with too high n proportion
	foreach my $loc(sort keys %ncount){
		if (($nAcount{$loc} / $admixinds) > $threshold){
			$$blacklistref{$loc} = ($nAcount{$loc} / $admixinds) unless exists $$blacklistref{$loc};
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
