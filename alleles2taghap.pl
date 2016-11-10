#! /usr/bin/perl

# Script by: Tyler K. Chafin
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
getopts( 'a:o:c:hp:x:i:m:sn:', \%opts );

# kill if help option is true
if( $opts{h} ){
  &help;
  die "Program killed: Help menu called.\n\n";
}

my $singleton = 0;
#Check boolean options
if( $opts{s} ){
  $singleton = 1;
}

#get options 
my ($in, $out, $cov, $popmap, $include, $exclude, $maxSNP, $indProp) = &parseArgs(\%opts); 

#Validate some options
$indProp > 1.0 and die "\nFatal error: -n cannot be greater than 1.0!\n\n";
$indProp < 0.0 and die "\nFatal error: -n cannot be less than 0.0!\n\n";
$maxSNP < 1 and die "\nFatal error: -m must be a positive integer!\n\n";
$cov < 0 and die "\nFatal error: -c must be an integer equal to or greater than 0!\n\n";

#Extract populations to include/exclude into array formats 
my (@popInc, @popExc) = ();
@popInc = split(/\+/,$include) unless !($include);
@popExc = split(/\+/,$exclude) unless !($exclude);

#Report user choices
print "\nParsing popmap...\n";
print "Populations to include: @popInc\n" if @popInc; 
print "Population to exclude: @popExc\n" if @popExc;
if (!@popInc && !@popExc){
	print "No populations specified, all samples in popmap will be used.\n";
}

#parse popmap file into a hash with ind as key and count as value
#only keep individuals from accepted pops
#Count (hash values) will be used to track number of loci found per sample
my $assignRef = &parsePopmapSelect($popmap, \@popInc, \@popExc); 
my $inds = keys %{$assignRef};
print "Total individuals retained from popmap: <$inds>\n";


print "\nParsing data file...\n";
$singleton == 1 and print "Skipping parsimony uninformative SNPs...\n";
print "Skipping loci with greater than <$maxSNP> SNPs...\n\n";

#Parse .alleles file, result is reference to array (each locus) of hashes (each sample) of arrays (each allele)
my $dataRef = &parseAllelesFile($in);

#Parse data structure, move loci with too low coverage, and individuals containing missing data or gaps
#Operates on array output from parseAllelesFile, passed by reference, returns similar struct but trimmed
#Function DOES alter $assignRef to add count data for each sample !!! 
my $cleanDataRef = &filterAllelesMatrix($dataRef, $assignRef, $cov, $singleton, $maxSNP);

#Report coverage per individual:
open (IND, ">$out.indcov") || die "Cannot open $out.indcov: $!\n";
print "\nReporting individual coverage in file $out.indcov...\n";
my $indCov = int(scalar @{$cleanDataRef} * $indProp);
my $indPercent = $indProp*100;
print "Retaining individuals with at least <$indPercent%> (N=<$indCov>) passing loci present...\n";
print IND "Individual	LociFound	Retained?\n";
my $kept = 0; 
for my $ind (sort keys $assignRef){
	print IND "$ind \t $assignRef->{$ind}\t";
	if ($assignRef->{$ind} >= $indCov){
		print IND "Y\n";
		$kept++; 
	}else{
		print IND "N\n";
		delete($assignRef->{$ind}); #Delete individuals having too low coverage
	}
}
close IND; 
print "Kept <$kept> individuals, out of a total of <$inds> included.\n";
$kept == 0 and die "\nProgram killed: Zero individuals kept after filtering!\n\n";

#Print output file, leaving empty for missing data
#And finally filter loci for coverage of individuals
print "\nRemoving loci with less than <$cov> individuals...\n";
#Pass datastructure (Array of loci, each locus is hash of individual alleles (Keys: Ind1_0, Ind1_1)
#AssignRef is hash of individual names which passed filtering
#$cov is number of individuals required to keep a locus
&printTagHap($cleanDataRef, $assignRef, $cov, $out);

print "\nFinal output saved in SimpleMatrix file <$out.taghap>!\n\n";

exit;

 ########################### SUBROUTINES ###############################

 sub help{
	 
	print "\nalleles2taghap.pl by Tyler Chafin\n";
	print "\nThis script converts from the .alleles file output by pyRAD to create the input for fineRADstructure\n";
	print "\nNOTE: 
	- All samples are assumed to be diploid.\n";
	print "	- Sample names CANNOT contain underscores.\n";
	print "	- Columns containing Ns or gaps will be deleted from final output\n";
	print "	- Popmap file should be tab-delimited, like so: SampleName [tab] PopID\n";
	print "	- If populations to include/exclude are not given, all samples in popmap are used.\n";
	print "	- You can specify multiple popIDs as: ID1+ID2+ID3, as long as these match IDs in popmap\n";
	print "	- For the -s filter, singletons are evaluated within the selected subset of individuals\n";
	print "	- Three-nucleotide ambiguities are treated as N's (V,H,D,B)\n";
	print "\nOptions:\n";
	print "\t-a	: Path to input file (.alleles)\n";
	print "\t-p	: Path to popmap file (tab-delimited)\n";
	print "\t-o	: Output file prefix. [Default = out, i.e. out.taghap]\n";
	print "\t-c	: Min number of samples for which data must be present per locus [Default = 1]\n";
	print "\t-n	: Minumum proportion of loci an individual must be present at to be retained [def = 0.2]\n";
	print "\t-i	: PopIDs to include in output file (e.g. -i pop1+pop4)\n";
	print "\t-x	: PopIDs to exclude (e.g. -x catenatus or -x sistrTX+sistrIN)\n";
	print "\t-m	: Maximum number of SNPs per locus. Loci exceeding are deleted [default:10]\n";
	print "\t-s	: Skip SNPs that are singletons [Boolean; Default = false]\n";
	print "\t-h	: Displays this help message\n";
	print "\n";
}

#parse arguments
sub parseArgs{

  my( $params ) =  @_;
  my %opts = %$params;
  
  #defaults
  my $in = $opts{a} or die "\nInput file not specified.\n\n";
  my $out = $opts{o} || "out"; 
  my $cov = $opts{c} || 1;
  my $map = $opts{p} or die "\nPopmap file not specified.\n\n";
  my $inc = $opts{i};
  my $exc = $opts{x};
  my $max = $opts{m} || 10;
  my $ind = $opts{n} || 0.2;
  #return
  return ($in, $out, $cov, $map, $inc, $exc, $max, $ind);
}

#Subroutine to parse .alleles file for variable columns only
sub parseAllelesFile{

  my $infile = $_[0];

  open(IN, $infile) || die "Failed to open input file: $!\n\n";
  
  my %temp;
  my @allLoci;
  my %locus; 
  my $linecount = 0;
  my $indcount = 0;
  my $loccount = 0;
  my $len;
  
  #For each line, collect data until end of locus, then parse locus (temp hash)
  while (<IN>){
    chomp; 
    if ($indcount == 0){
      $_ =~ /(\>.*\s+)/;
      my $match = $1;
      $match =~ s/ /@/g;
      #Get length of sample name, for use in parsing last line of locus
      $len = length($match);
      #print "Length: $len\n";
    }
    #Remove first character of name
    $_ =~ s/\>//;
    #split sample name and data, use to load temperary hash
    my @line = split(/\s+/,$_);
    
    #If this is end of locus, capture results
    if ($line[0] eq "//"){
      #Report num inds found this locus
      $loccount++;
      #print "Reading locus: $loccount, found <", ($indcount/2), "> individuals.\n";
      #Call function to get positions of SNPs
      my $snps = &parseLocusInfo($len, $_);
      #print "Line: $linecount, @{$snps}\n";
      #Use SNPs positions to subtr out the variable columns
      #Retain these in final hash
      foreach my $key (keys %temp){
        #print "$key:\n";
        foreach my $allele (@{$temp{$key}}){
          my $vars = "";
          foreach my $i (@{$snps}){
            $vars .= substr $allele, $i, 1;
          }
          #print "$vars\n";
          #Push string of SNPs ("haplotype") into hash
          push(@{$locus{$key}}, $vars);
        }
      }
      #Reset counters, move to next locus
      $indcount = 0;
      $len = 0;
      %temp = ();
      #This pushes an anonymous hash, a copy of current %locus hash, into allLoci arr
      push(@allLoci, {%locus});
      %locus = ();
      next;
    }
    my @name = split(/_/,$line[0]);
    #Push into temp array, to process when we have loaded whole locus
    push(@{$temp{$name[0]}}, $line[1]);
    #print "$name[0]: $temp{$name[0]}->[$#{$temp{$name[0]}}]\n";
    #Increment line count
    $indcount++;
    $linecount++;
  }

  close IN;
  print "Found <$loccount> loci.\n";
  return (\@allLoci);

}

#Subroutine for internal use, grabs SNP positions from
sub parseLocusInfo{
  my $offset = $_[0];
  my $string = $_[1];
  
  #To store positions
  my @temp; 
  my $pos;
  my $os = 0;
  #Get variable (not parsimony-informative) sites -- Comment out if these aren't wanted
  while (1){
    $pos = index ($string, "-", $os);
    last if ($pos < 0);
    push (@temp,  ($pos-$offset));
    $os = $pos+1;
  }
  
  #Get variable (parsimony-informative) sites
  $pos=-1;
  $os = 0;
    while (1){
    $pos = index ($string, "\*", $os);
    last if ($pos < 0);
    push (@temp,  ($pos-$offset));
    $os = $pos+1;
  }
  #Return sorted array
  return ([sort{$a <=> $b} @temp]);
}

#parse popmap file, filtering for selections
sub parsePopmapSelect{
	
	my $toParse = $_[0]; 
	my $inc = $_[1];
	my $exc = $_[2];
	
	#vars
	my %toReturn; 
	
	#Hash pop selections
	my %incMap = map { $_ => 0 } @{$inc};
	my %excMap = map { $_ => 0 } @{$exc};
	
	#open popmap
	open (POP, $toParse) or die "Oh no! Cannot open $toParse: $!\n"; 
	
	while (my $line = <POP>){
	  chomp $line; 

	  #ignore if blank
      if( $line =~ /^\w/ ){
        my @temp = split( /\s+/, $line); 
        
        #If includeMap is given, but pop is not included, skip
        if ((keys %incMap > 0) && !(exists $incMap{$temp[1]})){
          next;
        #else if exclude map given, and pop is excluded, skip
        }elsif ((keys %excMap > 0) && (exists $excMap{$temp[1]})){
          next; 
        #Otherwise, push ind into hash
        }else{
          #push into our hash
          $toReturn{$temp[0]} = 0;
        }
      }
	}
	
	close POP;
	return( \%toReturn);
}

#Subroutine to parse alleles hash, where indivs are hash keys, and each value is an 
#array of 2 "alleles" haplotypes, with just SNPs included.
sub filterAllelesMatrix{
  
  my $data 		= $_[0];
  my $counts	= $_[1];
  my $threshold	= $_[2];
  my $single 	= $_[3];
  my $max		= $_[4];

  #open(LOG, ">locus_log.tsv") || die "filterAllelesMatrix: Unable to open logfile <locus_log.tsv> for output- $!\n\n";
  
  #Track some stuff 
  my $locCount = 0; 
  my $keepCount = 0;
  my @finalMatrix; 
  
  #For each locus, grab alleles, check coverage, filter, and 
  foreach my $locus (@{$data}){
    my %tempMatrix;

    #Grab allele data for each individual for parsing
    foreach my $ind (sort keys %{$counts}){ 
      #If POPMAP ind has data for this locus, place in temp hash
      if (exists $locus->{$ind}){ 
        my $num = 0;
        #${$counts{$ind}}++; #Track how many loci an individual had data for
        foreach my $allele (@{$locus->{$ind}}){
           my $temp = ($ind."_".$num);
           #print $temp, " :\t";
           $tempMatrix{$temp} = $allele;
           #print "$allele\n";
           $num++;
        }
      }
    }
    $locCount++;
    #Check if locus has enough individuals to keep
    my $k = keys %tempMatrix;
    (($k/2) < $threshold) and next; #Skip if too few inds -> Thus not added to final data matrix
    #Get transposed matrix of allele sequences
    my $columnRef = getColumnsMatrix(\%tempMatrix);
    
    #Process locus to get columns to keep 
    my $keep = parseColumnMatrix($columnRef, $single, $max);
    
    #Loop through tempMatrix to extract data to keep, filter by cov etc.
    #print "Locus $locCount:\n";
    #If no columns pass, skip to next
    scalar @{$keep} == 0 and next; 
    #If too many SNPs, skip to next
    scalar @{$keep} > $max and next; 
    my %filteredTempMatrix; 
    foreach my $indAllele (sort keys %tempMatrix){
      my $passed = "";
      foreach my $column (@{$keep}){
		#for each slected column, substring the sequence. 
		$passed .= substr $tempMatrix{$indAllele}, $column, 1;
	  }
	  #print "$indAllele : $passed\n";
	  if ($indAllele =~ /_0/){
		my $basename = $indAllele;
		$basename =~ s/_0//g;
		$counts->{$basename}++; #Track which individuals have data
		#print "$basename has data\n";
	  }
	  $filteredTempMatrix{$indAllele} = $passed;
	}
    #Push into a structure to return to main program
	$keepCount++; #Increment if locus passed all filtering
	#This pushes an anonymous hash, a copy of current %locus hash, into allLoci arr
    push(@finalMatrix, {%filteredTempMatrix});
  }
  print "Retained <$keepCount> loci after filtering...\n";
  $keepCount == 0 and die "\nProgram killed: Zero loci kept after filtering!\n\n";
  #Close logfile
  #close LOG;
  return (\@finalMatrix);
}


#Internal subroutine to put sequence into format where array position is the 
#column in alignment, and value is a string of all data for that column, returns array
#More code borrowed/stolen from Steve
sub getColumnsMatrix{

  my $hashRef = $_[0];

  my @align;  
  #For each ind
  foreach my $key( sort keys %{ $hashRef } ){
    my $index = 0;
    my @seq = split( //, $hashRef->{$key}  );
    #for each nucleotide 
    foreach my $item( @seq ){
      $align[$index] .= $item;
      $index++;
    }
  }
  return( \@align );
}

#Internal subroutine to select out variable columns from the columnMatrix
sub parseColumnMatrix{

  my $colMat = $_[0];
  my $single = $_[1];
  my $max = $_[2];
  my @snps; 
  my $count = -1;

  foreach my $i (@{$colMat}){
    $count++;
    #If there are gaps, skip column  
    $i =~ tr/-// > 0 and next;

    #If there are Ns, skip column
    $i =~ tr/NnVvHhDdBb// > 0 and next;
    
    my $a = ($i =~ tr/Aa// > 0);
    my $c = ($i =~ tr/Cc// > 0);
    my $g = ($i =~ tr/Gg// > 0);
    my $t = ($i =~ tr/Tt// > 0);
    
    
    if ($a+$c+$t+$g <= 1){
      next; 
    }else{ 
      if ($single == 1){

        my @nums;
        #If evaluating singletons, get the number of freq of each allele
        $a == 1 and push (@nums, my $ai = ($i =~ tr/Aa//));
        $t == 1 and push (@nums, my $ti = ($i =~ tr/Tt//));
        $c == 1 and push (@nums, my $ci = ($i =~ tr/Cc//));
        $g == 1 and push (@nums, my $gi = ($i =~ tr/Gg//));
        my $alleles = scalar @nums;
        $alleles == 1 and next; #If only one allele, skip (monomorphic)
        #If there are two alleles, and one is a singleton, skip!
        if ($alleles == 2){
           $nums[0] == 1 and next; 
           $nums[1] == 1 and next; 
        }
       # print "Nums are @nums\n";
      }
      #Keep index of SNP
      push(@snps, $count);
    }
    #print "$i\n";
  }
  #print "@snps\n";
  return(\@snps);
}

sub printTagHap{
	
	my $data = $_[0];
	my $inds = $_[1];
	my $cov = $_[2];
	my $name = $_[3];
	
	my $kept = 0;
	
	#Open file for writing
	open (OUT, ">$name.taghap") || die "\nCannot open file $name.taghap: $!\n\n";
	
	#Print header line of individual labels
	#MUST retain this order for all other lines
	my $line = 0;
	for my $ind (sort keys %{$inds}){
		if ($line > 0){
			print OUT "\t";
		}
		print OUT "s", $ind;
		$line++; 
	}
	print OUT "\n"; 
	
	#Now print a line for each locus, leaving blank for missing data
	for my $loc (@{$data}){
		#Check coverage for locus
		my $temp = scalar(keys %{$loc}); #Number of alleles in hash
		if (($temp/2) < $cov){
			next; #Too few individuals present
		}else{
			$kept++; #Track number of kept loci
		}
		#Loop through individuals to fetch data
		my $count = 0;
		for my $ind (sort keys %{$inds}){
			my $a1 = $ind . "_0";
			my $a2 = $ind . "_1";
			$count > 0 and print OUT "\t";
			if (exists $loc->{$a1} && exists $loc->{$a2}){ 
				if ($loc->{$a1} eq $loc->{$a2}){
					print OUT $loc->{$a1};
				}else{
					print OUT $loc->{$a1}, "/", $loc->{$a2};
				}
			}
			$count++; 
		}
		print OUT "\n";
	}
	close OUT;
	$kept == 0 and die "\nProgram killed: Zero loci kept after checking for individual coverage!\n\n";
	print "Retained <$kept> loci after filtering for coverage...\n";
}













