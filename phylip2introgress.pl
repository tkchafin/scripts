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
getopts( 'p:i:1:2:a:o:hn:N:gPxrOd:', \%opts );

# kill if help option is true
if( $opts{h} ){
  &help;
  die "Help menu\n\n";
}

#get options
my ($map, $phy, $p1, $p2, $ad, $out, $threshold, $globalThresh, $gapFalse, $phyNew, $onlyPhy, $rcalls, $order, $delta) = &parseArgs(\%opts);

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
#gapFalse = 1 only count Ns; otherwise count gaps
&getBlacklist($threshold, $globalThresh, $pop1Ref, $pop2Ref, $popaRef, $nchar, $gapFalse, \%blacklist);

my $nFail = (keys %blacklist);

print($nFail ," loci had greater than ",$threshold, " missing data. Removing them.\n");
#print Dumper(\%blacklist);

if ($delta > 0.0){
	&getDeltaBlacklist($delta, $pop1Ref, $pop2Ref, $nchar, \%blacklist);
	my $nDelt = (keys %blacklist);
	my $diff = $nDelt - $nFail;
	print($diff ," loci had delta less than ",$delta, ". Removing them.\n");
}

$nFail = (keys %blacklist);
if ($nFail >= $nchar){
	print("Oh no! All loci failed filtering! Try different values of -n, -N, or -d!\n");
	exit;
}else{
	my $nPass = $nchar - $nFail;
	print($nPass, " loci passed filtering!\n")
}

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
if ($onlyPhy != 1){
	open(ADMIX, "> admix.csv");
	open(LOCI, "> loci.txt");
	open(P1DATA, "> p1data.csv");
	open(P2DATA, "> p2data.csv");

	#Print loci.txt file
	print LOCI "locus,type\n";
	for my $nloci (0 .. $nchar-1){
		if(!exists $blacklist{$nloci}){
			print LOCI "loci$nloci,C\n";
		}
	}

	print "Done writing LOCI file <loci.txt>\n";
	close LOCI;

	#Print header lines of admix file
	my @line1;
	my @line2;
	foreach my $adInd (keys %{$popaRef}){
		push( @line1, "pop1");
		push( @line2, $adInd);
	}

	if ($order == 1){
		open (ORDER, ">pop_order_admix.txt");
		foreach my $adInd (sort keys %{$popaRef}){
			print ORDER $adInd, "\t", ($assignRef->{$adInd}, "\n");
		}
		close(ORDER)
	}

	my $line1str = join(",", @line1);
	my $line2str = join(",", @line2);

	print ADMIX $line1str, "\n";
	print ADMIX $line2str, "\n";

	#format and print files for introgress
	for (my $loc = 0; $loc < $nchar; $loc++){
		if(!exists $blacklist{$loc}){
			#write p1 file
			my @p1line;
			foreach my $ind (sort keys %{$pop1Ref}){
				my $nuc = ${$pop1Ref->{$ind}}->[$loc];
				$nuc =~ s/A/A\/A/g;
				$nuc =~ s/T/T\/T/g;
				$nuc =~ s/G/G\/G/g;
				$nuc =~ s/C/C\/C/g;
				$nuc =~ s/W/A\/T/g;
				$nuc =~ s/R/A\/G/g;
				$nuc =~ s/M/A\/C/g;
				$nuc =~ s/K/G\/T/g;
				$nuc =~ s/Y/T\/C/g;
				$nuc =~ s/S/C\/G/g;
				$nuc =~ s/N/NA\/NA/g;
				$nuc =~ s/-/NA\/NA/g;
				push(@p1line, $nuc);
			}
			my $p1linestr = join(",", @p1line);
			print P1DATA $p1linestr, "\n";

			#Write p2 file
			my @p2line;
			foreach my $ind (sort keys %{$pop2Ref}){
				my $nuc = ${$pop2Ref->{$ind}}->[$loc];
				$nuc =~ s/A/A\/A/g;
				$nuc =~ s/T/T\/T/g;
				$nuc =~ s/G/G\/G/g;
				$nuc =~ s/C/C\/C/g;
				$nuc =~ s/W/A\/T/g;
				$nuc =~ s/R/A\/G/g;
				$nuc =~ s/M/A\/C/g;
				$nuc =~ s/K/G\/T/g;
				$nuc =~ s/Y/T\/C/g;
				$nuc =~ s/S/C\/G/g;
				$nuc =~ s/N/NA\/NA/g;
				$nuc =~ s/-/NA\/NA/g;
				push(@p2line, $nuc);
			}
			my $p2linestr = join(",", @p2line);
			print P2DATA $p2linestr, "\n";

			#Populate admix file
			my @admixline;
			foreach my $ind (sort keys %{$popaRef}){
				my $nuc = ${$popaRef->{$ind}}->[$loc];
				$nuc =~ s/A/A\/A/g;
				$nuc =~ s/T/T\/T/g;
				$nuc =~ s/G/G\/G/g;
				$nuc =~ s/C/C\/C/g;
				$nuc =~ s/W/A\/T/g;
				$nuc =~ s/R/A\/G/g;
				$nuc =~ s/M/A\/C/g;
				$nuc =~ s/K/G\/T/g;
				$nuc =~ s/Y/T\/C/g;
				$nuc =~ s/S/C\/G/g;
				$nuc =~ s/N/NA\/NA/g;
				$nuc =~ s/-/NA\/NA/g;
				push( @admixline, $nuc);
			}
			my $admixlinestr = join(",", @admixline);
			print ADMIX $admixlinestr, "\n";
		}
	}

	print "Done writing P1DATA file <p1data.csv>\n";
	close P1DATA;
	print "Done writing P2DATA file <p2data.csv>\n";
	close P2DATA;
	print "Done writing ADMIX file <admix.csv>\n\n";
	close ADMIX;
}

if ($phyNew or $onlyPhy){
	print("Writing new PHYLIP file <out.phy>\n");
	open (PHY, "> out.phy");
	my $locnum = 0;
	my $indnum = 0;
	for (my $loc = 0; $loc < $nchar; $loc++){
		if(!exists $blacklist{$loc}){
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
			if(!exists $blacklist{$l}){
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
			if(!exists $blacklist{$loc}){
				print PHY ${$pop2Ref->{$ind}}->[$loc];
			}
		}
		print PHY "\n";
	}
	foreach my $ind (sort keys %{$popaRef}){
		print PHY $ind, "\t";
		for (my $loc = 0; $loc < $nchar; $loc++){
			if(!exists $blacklist{$loc}){
				print PHY ${$popaRef->{$ind}}->[$loc];
			}
		}
		print PHY "\n";
	}
	close PHY;
}

#if $rcall, print some stuff
if ($rcalls == 1){
	print"Printing some population-wise R calls....\n";

	my %samples;

	foreach my $ind (keys %{$pop1Ref}){
		if (!exists $samples{$assignRef->{$ind}}){
			$samples{$assignRef->{$ind}} = [$ind];
		}else{
			push(@{$samples{$assignRef->{$ind}}}, $ind);
		}
	}
	foreach my $ind (keys %{$pop2Ref}){
		if (!exists $samples{$assignRef->{$ind}}){
			$samples{$assignRef->{$ind}} = [$ind];
		}else{
			push(@{$samples{$assignRef->{$ind}}}, $ind);
		}
	}
	foreach my $ind (keys %{$popaRef}){
		if (!exists $samples{$assignRef->{$ind}}){
			$samples{$assignRef->{$ind}} = [$ind];
		}else{
			push(@{$samples{$assignRef->{$ind}}}, $ind);
		}
	}

	foreach my $pop (keys %samples){
		print "$pop<- genomic.clines(introgress.data=count.matrix, hi.index=hi.index.sim,loci.data=loci.data.sim,sig.test=T,classification=T, loci.touse=which(deltas>0.8), ind.touse=c(";
		my $count = 0;
		foreach my $assign (@{$samples{$pop}}){
			if ($count == 0){
				print "\"$assign\"";
			}else{
				print ", \"$assign\"";
			}
			$count++;
		}
		print "))\n";
	}
	print "\n";
}



exit 0;

 ########################### SUBROUTINES ###############################

 sub help{

	print "\nphylip2introgress.pl by Max Bangs, Tyler Chafin and Steve Mussmann\n";
	print "\nThis script converts from phylip format to the input for Introgress R package.\n";
	print "A population map should be given in a tab-delimited file, formatted as:\n";
	print "\n\tSampleName\tPopID\n\n";
	print "Where PopID can be a string or integer, and SampleName must exactly match a corresponsing line in the phylip file. Any samples not in the popmap will not be included in the output files.\n\n";
	print "Options:\n";
	print "\t-p	: Path to popmap file (tab-delimited)\n";
	print "\t-1	: Identifier for population 1 (include multiple as: pop1+pop2)\n";
	print "\t-2	: Identifier for population 2 (include multiple as: pop1+pop2)\n";
	print "\t-a	: Identifier for admixed population(s) (include multiple as: pop1+pop2)\n";
	print "\t-i	: Path to input file (phylip)\n";
	print "\t-n	: Proportion missing data allowed per population per SNP (default=0.5)\n";
	print "\t-N	: Proportion of globally missing data allowed per SNP (default=0.5)\n";
	print "\t-g	: Toggle on to TURN OFF default behavior of treating gaps as missing data\n";
	print "\t-P	: Toggle on to output a new phylip file with the filtered data. [default=off]\n";
	print "\t-x	: Toggle on to ONLY print a new phylip file (e.g. no INTROGRESS files)\n";
	print "\t-r	: Toggle on to print a file with population-wise genomic.clines calls\n";
	print "\t-d	: Delta threshold to retain a locus [Default=0.0]\n";
	print "\t(Delta = sigma(|fi1 - fi2|/2) where fi1=freq of ith allele in pop 1; Gregorius & Roberds 1986)\n";
	print "\t-O	: Toggle on to print a file of pop assignments for ADMIX pops <pop_order_admix.txt>\n";
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
	$opts{x} and $onlyPhy = 1;
  my $out = $opts{o} || "out.phy";
	my $r = 0;
	$opts{r} and $r = 1;
	my $ord = 0;
	$opts{O} and $ord = 1;
	my $del = $opts{d} || 0.0;
  #return
  return ($map, $phy, $p1, $p2, $ad, $out, $threshold, $globalThresh, $gapFalse, $phyNew, $onlyPhy, $r, $ord, $del);
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
        my @arr = split(//, uc($temp[1]));
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

sub getDeltaBlacklist{
	my ( $delta_thresh, $p1ref, $p2ref, $locnum, $blacklistref ) = @_;
	print "Delta threshold is: $delta_thresh\n";
	for my $loc (0 .. $nchar-1){
		#skip if locus already in blacklist
		
		if(!exists $blacklistref->{$loc}){
			#print("Locus $loc\n");
			my %p1counts = (
				"A" => 0,
				"T" => 0,
				"G" => 0,
				"C" => 0
			);
			my %p2counts = (
				"A" => 0,
				"T" => 0,
				"G" => 0,
				"C" => 0
			);
			my $p1inds = 0;
			#get allele counts for p1 and p2 pops
			foreach my $ind( sort keys %$p1ref ){
				my $nuc = ${$$p1ref{$ind}}->[$loc];
				#print("$nuc\n");
				if ($nuc eq "A") { $p1counts{"A"} += 2;$p1inds++;} #allele count ONLY counts acceptable diploid genotypes
				elsif ($nuc eq "T") { $p1counts{"T"} += 2;$p1inds++;}
				elsif ($nuc eq  "C") { $p1counts{"C"} += 2;$p1inds++;}
				elsif ($nuc eq  "G") { $p1counts{"G"} += 2;$p1inds++;}
				elsif ($nuc eq  "R") { $p1counts{"A"}++; $p1counts{"G"}++;$p1inds++;}
				elsif ($nuc eq  "Y") { $p1counts{"C"}++; $p1counts{"T"}++;$p1inds++;}
				elsif ($nuc eq  "S") { $p1counts{"G"}++; $p1counts{"C"}++;$p1inds++;}
				elsif ($nuc eq  "W") { $p1counts{"A"}++; $p1counts{"T"}++;$p1inds++;}
				elsif ($nuc eq  "K") { $p1counts{"G"}++; $p1counts{"T"}++;$p1inds++;}
				elsif ($nuc eq  "M") { $p1counts{"C"}++; $p1counts{"A"}++;$p1inds++;}
			}
			my $p2inds=0;
			foreach my $ind( sort keys %$p2ref ){
				my $nuc = ${$$p2ref{$ind}}->[$loc];
				#print("$nuc\n");
				if ($nuc eq "A") { $p2counts{"A"} += 2;$p2inds++;}
				elsif ($nuc eq "T") { $p2counts{"T"} += 2;$p2inds++;}
				elsif ($nuc eq  "C") { $p2counts{"C"} += 2;$p2inds++;}
				elsif ($nuc eq  "G") { $p2counts{"G"} += 2;$p2inds++;}
				elsif ($nuc eq  "R") { $p2counts{"A"}++; $p2counts{"G"}++;$p2inds++;}
				elsif ($nuc eq  "Y") { $p2counts{"C"}++; $p2counts{"T"}++;$p2inds++;}
				elsif ($nuc eq  "S") { $p2counts{"G"}++; $p2counts{"C"}++;$p2inds++;}
				elsif ($nuc eq  "W") { $p2counts{"A"}++; $p2counts{"T"}++;$p2inds++;}
				elsif ($nuc eq  "K") { $p2counts{"G"}++; $p2counts{"T"}++;$p2inds++;}
				elsif ($nuc eq "M") { $p2counts{"C"}++; $p2counts{"A"}++;$p2inds++;}
			}
			#exit;
			my $delta = 0.0;
			foreach my $key( sort keys %p1counts ){
				my $p1_freq = $p1counts{$key}/($p1inds*2); #allele count / 2N
				my $p2_freq = $p2counts{$key}/($p2inds*2);
				my $diff = abs($p1_freq - $p2_freq) / 2.0;

				$delta += $diff;
			}
			if ($delta <= $delta_thresh){
				#print("Delta:$delta\n");
				$$blacklistref{$loc} = "delta" unless exists $$blacklistref{$loc};
				#print("Locus:",$loc," - Missing data: ",($n2count{$loc} / $p2inds), "\n");
			}
		}
	}
	#exit;
}

sub getBlacklist{

	my( $thresh, $globalThresh, $p1ref, $p2ref, $admixref, $nchar, $gap, $blacklistref ) = @_;

	my $globalInds = 0;
	my %globalCount;

	#Check loci in pop1
	my %ncount;
	my $p1inds;
	foreach my $ind( sort keys %$p1ref ){
		$p1inds++;
		$globalInds++;
		for my $loc (0 .. $nchar-1){
			my $nuc = ${$$p1ref{$ind}}->[$loc];
			$ncount{$loc} = 0 unless exists $ncount{$loc};
			$globalCount{$loc} = 0 unless exists $globalCount{$loc};
			if($gap==1 and $nuc eq "N"){
				$ncount{$loc}++;
				$globalCount{$loc}++;
			}
			if($gap==0){
				if($nuc eq "N" or $nuc eq "-"){
					$ncount{$loc}++;
					$globalCount{$loc}++;
				}
			}
		}
	}
	#blacklist any loci with too high n proportion
	foreach my $loc(sort keys %ncount){
		if (($ncount{$loc} / $p1inds) > $threshold){
			#print("$loc: $ncount{$loc}\n");
			foreach my $ind( sort keys %$p1ref ){
				my $nuc = ${$$p1ref{$ind}}->[$loc];
				#print("$nuc\n");
			}
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
		for my $loc (0 .. $nchar-1){
			my $nuc = ${$$p2ref{$ind}}->[$loc];
			$n2count{$loc} = 0 unless exists $n2count{$loc};
			$globalCount{$loc} = 0 unless exists $globalCount{$loc};
			if($gap==1 and $nuc eq "N"){
				$n2count{$loc}++;
				$globalCount{$loc}++;
			}
			if($gap==0){
				if($nuc eq "N" or $nuc eq "-"){
					$n2count{$loc}++;
					$globalCount{$loc}++;
				}
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
		#print("$ind\n");
		for my $loc (0 .. $nchar-1){
			my $nuc = ${$$admixref{$ind}}->[$loc];
			#print("$nuc\n");
			$nAcount{$loc} = 0 unless exists $nAcount{$loc};
			$globalCount{$loc} = 0 unless exists $globalCount{$loc};
			if($gap==1 and $nuc eq "N"){
				$nAcount{$loc}++;
				$globalCount{$loc}++;
			}
			if($gap==0){
				if($nuc eq "N" or $nuc eq "-"){
					$nAcount{$loc}++;
					$globalCount{$loc}++;
				}
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
