#! /usr/bin/perl

#Written by Steven M. Mussmann 
#TKC Modifications:
#17-July-2016 
#  -r: Random sample SNPs
#December 2015
#  -s: Filter singleton SNP loci 
#  -b: Retain only biallelic SNPs (needed for SNAPP and other analyses) 
#  - Minor bug fix for MAF filter (option a) 
#  - "filter" subroutine: Consider biallelic ambiguities valid, whcih assumes diploid alignment

use warnings;
use strict;
use Getopt::Std;
#use Data::Dumper;

# kill program and print help if no command line arguments were given
if( scalar( @ARGV ) == 0 ){
  &help;
  die "Exiting program because no command line options were used.\n\n";
}

# take command line arguments
my %opts;
getopts( 'a:bc:f:g:hi:mo:st:r:', \%opts );

# if -h flag is used, or if no command line arguments were specified, kill program and print help
if( $opts{h} ){
  &help;
  die "Exiting program because help flag was used.\n\n";
}
 
#lol 
&stupidpic;

# parse the command line
my( $file, $individual, $column, $out, $type, $minfreq, $remMonomorph, $remSingleton, $biSNP, $gaps, $rand ) = &parsecom( \%opts );

# declare global variables
my $hashref; # hashreference that will hold sample and sequence data that is returned from subroutines that parse input files

# parse the file
if( $type =~ /fasta/i ){
  $hashref = &fastaparse( $file );
}elsif( $type =~ /phylip/i ){
  $hashref = &phyparse( $file );
}else{
  die "File type $type is unrecognized.  Only Fasta and Phylip are supported.\n\n";
}

# find list of keys to be removed
my( $keyremove ) = &filter( $hashref, 0, $individual, 0, 0, 1, 0, 1.0);

# delete any sequence from the hash that is included in the $keyremove list
foreach my $item( @$keyremove ){
  delete $$hashref{$item};
}

# copy data into $alignref with the index value of the alignment as the key and a string of nucleotides at that index as the value
my( $alignref ) = &getcolumns( $hashref );

#Parse for columns to remove due to not meeting thresholds 

my $columnremove = &filter( $alignref, $minfreq, $column, $remMonomorph, $remSingleton, 0, $biSNP, $gaps);

#remove columns flagged as too low MAF or too high ambiguity
&removecolumns( $hashref, $columnremove ); 

#If random sample invoked
if ( $rand ){
  my ( $tempref ) = &getcolumns( $hashref );
  my ( $cull ) = &randomSample( $tempref, $rand );
  &removecolumns ( $hashref, $cull );
}

# print the output
if( $type =~ /fasta/i ){
  &fastaprint( $out, $hashref );
}elsif( $type =~ /phylip/i ){
  &phyprint( $out, $hashref );
}else{
  die "File type $type is unrecognized.  Only Fasta and Phylip are supported.\n\n";
}

# print success message to STDOUT
print "\nOutput written to $out\n\n";

#Report missing data content
&reportMatrixContent($hashref);

#print Dumper( \%$hashref );
#print Dumper( \%$alignref );

exit;

#####################################################################################################
############################################ Subroutines ############################################
#####################################################################################################

# subroutine to print help
sub help{
  
  print "\nnremover.pl is a perl script developed by Steven Michael Mussmann\n\n";
  print "Last Modified: Tyler K. Chafin - 17-July-2016; Added random sample option\n\n";
  print "To report bugs send an email to mussmann\@email.uark.edu\n";
  print "When submitting bugs please include all input files, options used for the program, and all error messages that were printed to the screen\n\n";
  print "Program Options:\n";
  print "\t\t[ -a | -b | -c | -f | -h | -i | -m | -o | -s | -t | -r ]\n\n";
  print "\t-c:\tUse this flag to specify the maximum proportion of ambiguities allowed in a column of a DNA sequence alignment.\n";
  print "\t\tDefault value is 0.2.\n\n";
  print "\t-f:\tUse this flag to specify the input file name.\n";
  print "\t\tThere is no default value.  The program will terminate if no file is provided.\n\n";
  print "\t-h:\tUse this flag to display this help message.\n";
  print "\t\tThe program will die after the help message is displayed.\n\n";
  print "\t-i:\tUse this flag to specify the max proportion of ambiguities allowed in an individual.\n";
  print "\t\tDefault value is 0.5.\n\n";
  print "\t-o:\tUse this flag to specify the output file name.\n";
  print "\t\tIf no name is given, the file extension \".out\" will be appended to the input filename.\n\n";
  print "\t-t:\tUse this flag to specify the input file type.\n";
  print "\t\tThe defualt file type is fasta.  Only fasta and phylip are currently supported.\n\n";

  print"\t-m:\tBoolean. Use this flag to toggle on removal of monomorphic loci.\n"; 
  print"\t\tThe default is false\n\n"; 

  print"\t-s:\tBoolean. Use this flag to toggle on removal of sites where the only SNP is a singleton.\n"; 
  print"\t\tThe default is false\n\n";

  print"\t-b:\tBoolean. Use this flag to only retain biallelic SNPs.\n";
  print"\t\tThe default is false. Using this flag forces [-m] option.\n\n";

  print"\t-a:\tUse this flag to specify minimum minor allele frequency.\n"; 
  print "\t\tThe default is 0.0. If any value is given, monomorphic loci will be deleted\n"; 
  print "\t\tGaps are treated as alleles, missing data are ignored.\n\n"; 

  print "\t-g:\tUse this flag to specify a maximum proportion of gaps allowed in a column.\n";
  print "\t\tThe default value is 1.0, or that any proportion of gaps is acceptable.\n\n"; 
  
  print "\t-r:\tUse this flag to subsample snps from the alignment (post-filtering).\n";
  print "\t\tThe default is to not subsample. Provide an integer value.\n\n"; 
  
  
}

#####################################################################################################
# subroutine to parse the command line options

sub parsecom{ 
  
  my( $params ) =  @_;
  my %opts = %$params;
  
  # set default values for command line arguments
  my $file = $opts{f} or die "\nNo input file specified\n\n"; # used to specify input file name.  Program will die if no file name is provided.
  my $individual = $opts{i} || 0.5; # used to specify threshold of ambiguities (Ns and other ambiguity codes) for each individual.  By default, if more than 50% of nucleotides are ambiguities for an individual, it will be discarded.
  my $column = $opts{c} || 0.2; # used to specify threshold of ambiguities (Ns and other ambiguity codes) for each column of an alignment.  By default, if more than 20% of nucleotides in a column are ambiguities, it will be discarded.
  my $out = $opts{o} || "$file.out"; #used to specify output file name.  If no name is provided, the file extension ".out" will be appended to the input file name.
  my $type = $opts{t} || "fasta"; #used to specify input file type.  Fasta is the default file type.  Only Fasta and Phylip are currently supported.
  
  my $minfreq = $opts{a} || 0.0; #Minimum minor allele frequency by default set to 0; meaning even monomorphic loci will be kept 

  my $remM = 0;
  my $remS = 0;  
  my $remB = 0;
  $opts{m} and $remM = 1; 
  $opts{s} and $remS = 1;  
  $opts{b} and $remB = 1;

  my $gapfreq = $opts{g} || 1.0;

  if ($opts{r}){ 
    !($opts{r} > 0) and die "\nInvalid use of [-r]!\n\n";
  }
  
  my $rand = $opts{r} || undef; #SNPs to random sample 
  
  return( $file, $individual, $column, $out, $type, $minfreq, $remM, $remS, $remB, $gapfreq, $rand );

}

#####################################################################################################
# subroutine to parse a fasta file

sub fastaparse{
  my( $file ) = @_;
  
  # declare strings to hold header and sequence
  my @head; # holds fasta headers
  my @seq; # holds sequences
  my %hash; # associates fasta headers with sequence data

  # open the fasta file
  open( FASTA, $file ) or die "Can't open $file, d-bag: $!\n";

  # loop through fasta file and extract data
  while( my $line = <FASTA> ){
    chomp( $line );
    if( $line =~ /^>/ ){
      push( @head, $line );
    }elsif( $line =~ /^[a-z]/i ){
      push( @seq, $line );
    }
  }
  
  for( my $i=0; $i<@seq; $i++ ){
    $hash{ $head[$i] } = uc($seq[$i]);
  }
    
  # delete original arrays
  undef @head;
  undef @seq;
  
  close FASTA;

  return( \%hash );

}

#####################################################################################################
# subroutine to parse a phylip file

sub phyparse{
  my( $file ) = @_;
  
  # declare strings to hold header and sequence
  my @lines;
  my @names;
  my @seq;
  my %hash;
  
  # open the phylip file
  open( PHYLIP, $file ) or die "Can't open $file, d-bag: $!\n";
  
  # loop through phylip file and put lines into array
  while( my $line = <PHYLIP> ){
    chomp( $line );
    # ignore blank lines
    if( $line =~ /^\w/ ){
      push( @lines, $line );
    }
  }

  # close file
  close PHYLIP;
  
  # Fix all whitespace in file so it is represented by a single tab
  for( my $i=0; $i<@lines; $i++ ){
    $lines[$i] =~ s/\s+/\t/g;
    # print $lines[$i], "\n";
  }
  
  # split each line to separate name and sequence into separate files
  foreach my $item( @lines ){
    my @temp = split( /\t/, $item );
    push( @names, $temp[0] );
    push( @seq, $temp[1] );
  }
  
  # remove first element from each array (this represents number of samples and nucleotides)
  shift( @names );
  shift( @seq );

  # put data into hash
  for( my $i=0; $i<@seq; $i++ ){
    $hash{ $names[$i] } = uc($seq[$i]);
  }
    
  # delete original arrays
  undef @names;
  undef @seq;
    
  # uncomment line to make sure hash is being created correctly
  #print Dumper( \%hash );

  return( \%hash );

}


####################################################################################################
# subroutine to remove columns with minor allele frequency below some value 
sub filter{ 

  my ($input, $threshold, $percent, $mono, $single, $ind, $bi, $g) = @_; 
  my @remove; #Cols to remove 

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
  
  #If freq too low, blacklist hash key 
  foreach my $key( sort keys %{$input} ){ 
    my $length = length(${$input}{$key});  
    if ($ind == 0){
 
      my %temp; #Store snp frequencies 
      my @snps = split //, ${$input}{$key}; #To iterate through the snps 

      foreach (@snps){ 
        if (my $tempkey = $codeHash{uc($_)}){
          foreach my $nuc(@{$tempkey}){ 
            $temp{$nuc}++; 
          } 
        } 
      }
      #Blacklist monomorphic loci 
      #Ambiguity >3 possible nucleotides are ignored in this 
      #Gaps count as an allele 
      my $size = keys %temp; 
      if ($mono == 1){
        if ($size <= 1){ 
          push @remove, $key; 
          next; 
        }
      }
      
      #Find key with smallest value; this is the minor allele
      if (my $minKey = minValue(\%temp)){   
        #If remove singleton option toggled
        if ($single == 1){
          #If there are exactly two alleles
          if ($size == 2){
            #And minor allele has frequency of 1, than it is singleton = remove 
            if ($temp{$minKey} == 1){
              push @remove, $key; 
              next; 
            }
          }else{
            #If num alleles != 2, and biallelic filter turned ON, remove locus
            if ($bi == 1){ 
              push @remove, $key; 
              next; 
            }
          }
        }

        #Calculate allele frequency of minor allele
        my $alleleFreq = ($temp{$minKey} / (2*$length));
        if (sprintf ("%.2f" ,$threshold) > 0.0){
          if (sprintf( "%.2f",  $alleleFreq) < sprintf("%.2f", $threshold)){ 
            push @remove, $key;   
            next; 
          }
        }
      }
    }

    #Filter by N content
    my $count = (${$input}{$key} =~ tr/NVHDB/NVHDB/ );
    my $ncontent = ( $count / $length );
    # fuck floating-point numbers
    if( sprintf( "%.2f", $ncontent ) > sprintf( "%.2f", $percent ) ){
       push @remove, $key;
       next; 
    } 

    #Filter for gaps content 
    if (sprintf("%.2f", $g) < 1.0){
      $count = (${$input}{$key} =~ tr/-/-/ ) || 0;
      my $gcontent = ( $count/ $length);
      if (sprintf( "%.2f", $gcontent) > sprintf( "%.2f", $g)){
        push @remove, $key;
        next; 
      }    
    }           
  }
  return (\@remove); 
}

sub minValue(){ 
  my $hash = shift; 
  keys %$hash; 
 
  my($minKey, $minVal) = each %$hash; 

  while (my ($key, $val) = each %$hash){ 
    if ($val < $minVal){ 
      $minVal = $val; 
      $minKey = $key; 
    } 
  }
  
  return $minKey; 
}


#####################################################################################################
# subroutine to print data out to a fasta file

sub fastaprint{

  my( $out, $hashref ) = @_;

  # open the output file for printing
  open( OUT, '>', $out ) or die "Can't open $out, d-bag: $!\n\n";

  # print the hash
  foreach my $key( sort keys %{ $hashref } ){
    print OUT $key, "\n";
    print OUT ${$hashref}{$key}, "\n";
  }

  # close the output file
  close OUT;

}

#####################################################################################################
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

#####################################################################################################
# subroutine to put sequence alignment into a hash with the index value of the alignment as the key and a string of nucleotides at that index as the value

sub getcolumns{

  my( $hashref ) = @_;
  
  my %align; # hash of arrays to hold position in alignment (hash key) and all characters at that position in alignment (value)
  
  
  foreach my $key( sort keys %{ $hashref } ){
    my $index = 0;
    my @sequence = split( //, ${$hashref}{$key}  );
    foreach my $item( @sequence ){
      $align{$index} .= $item;
      $index++;
    }
  }
  
  return( \%align );

}

#####################################################################################################
# subroutine to remove columns from an alignment, given the alignment contained in a hash and an array of positions in each value to be removed

sub removecolumns{

  my( $hashref, $columnremove ) = @_;

  # replace columns to be removed with a special character
  foreach my $key( sort keys %{ $hashref } ){
    foreach my $item( @$columnremove ){
      substr(${$hashref}{$key}, $item, 1) = 'z';
    }
  }
  
  # replace the special characters with nothing
  foreach my $key( sort keys %{ $hashref } ){
    ${$hashref}{$key} =~ s/z//g;
  }

}

#####################################################################################################
# subroutine to randomly sample SNPs to keep in alignment 

sub randomSample{

  my( $hashref, $num ) = @_;
  my @remove;
  
  my $total = scalar keys %$hashref;
  
  if ($num >= $total){
	print "\n[-r] invoked, but number to sample is greater than total SNPs remaining after filtering!\n\n";
    return (\@remove);
  }
  
  my $tokill = $total - $num;
 
  #Build array of indices
  my @indices = (0..$total-1);
  shuffle(\@indices); #Shuffle indices
  @remove = @indices[0..$tokill-1]; #Choose n of them for removal

  return (\@remove); #Return array of elements to remove

}

sub shuffle{
  my $array = shift;
  my $i = @$array; 
  while (--$i){
    my $j = int rand($i+1);
    @$array[$i,$j] = @$array[$j,$i];
  }
}

#####################################################################################################
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

#####################################################################################################
# prints a stupid ascii image just to be annoying
sub stupidpic{
	print '                                                                                        
                                                      T                                            
                                                     M                                             
                                              +MI   :                                              
                                                  ?.M ?MM888DMD,                                   
                                                   M88888888888888DN.                              
                                                 M88888888888888888888N~                           
                                                8888888888888888888888888M                         
                                              ,88888888888888888888888888888O                      
                                             .88888888888888888888888888888888M                    
                                            =88888888888888888888888888888888888M                  
                                           M88888888888888888888888888888888888888M                
                                         D888888888N88888888888888888888888888888888O              
                                             M8888MID888888888888888888888888N8888888D             
                                            M8888MIIM8888M8888888888888888888MIIMD8888D.           
                                           7888MDIIIM888MII888888888888888888IIIIIN ~M8M           
                                           = IIIIIID88TTIIIM888DI8888M888888?IIIIIII.              
                                             MIIIMIIMIIIII?M88M ,888IIMDIDDIIIM      ?             
                                             M?IIIIIIIIII.      M8MMIIIIIIIII.        M            
                                             MINIIIIIIII       ,    NIIIIIII$     R    N           
                                          .MI?N8IIIIIIIN      &7    7IIIIIIM      M7  M           
                                         .IIIIIIIMIIIII=      &&     DIIIIIIIIM        D           
                                         MIIIIIIIIIIIIIM             ?IIIIIIIIIIM     M            
                                         ?IIIDIIIIIIIIIIN           MIIIIIIIIIIIIIN M.             
                                         IIIIIOIIIIIIIIIIM        ,IIIIIIIIIIIIIIIII$              
                                         7IIIIIMIIIIIIIIIII?NMMM8IMIIIIIIIIIIIIIIIIIIM             
                                         MIIIIIIIIIIIIIIIIMMMMM7IIIIIIIIIIIIIIIIIIIIIID            
                                         .IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIM            
                                          NIIIIIIIIIIIIIIMIIIIIIIIIIIIIIIIIIIIMM?IIMM             
                                           NIIIIIIIIIIII7IIIIIIIIIMIIIIIIIIIIIIIIIIIIID.           
                                             IMMNIIIIIIIIIIIIIIII7IIIIMIIIMIIIIIII?II?IIIIIII8.    
                                                IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOIII?IIIIIIIIIII    
                                                IIIIIIIIINCLETUSD7IIIIIIIIIIIIIIIIIIIIIIIIIIIIM    
                                                IIIIIIMMMMMMM      O   IMM7IIIIIIIIIIIIIIIIID      
                                                IIIIMMMMMMMMMMM::MIM     N   .OMMZIIIIIDN,         
                                               :IIIMMMMMMMMMMMMIIIIM.,      M     ..             
                                               NIIINZDOZZMMMMMDIII~             MN~               
                                               IIIIMZZZZZZMMMIIIM                                 
                                              DIIIIIIMMMMZIIIIIIIIM                                
                                              IIIIIIIIIIIIIIIIIIIII,                               
                                             OIIIIIIIIIIIIIIIIIIIIN                                
                                             IIIIIIIIIIIIIIIIID..                                  
                                            MIIIIIIIIIIIIIIII,                                     
                                          ,IIIIIIIIIIIIIIIII?                                      
                                        DIIMIIIIIIIIIIIIIIIN                                       
                                    DOZZ7IIMIIIIIIIIIIIIIIO                                        
                                :M7MZZZLOLI7IIIIIIIIIIIIII~                                        
                             OYIIIIIZZZZMIIIIZIIIIIIIIIIIMN                                        
                         ,MIIIIIIIIIMZZZZMIIIIIMIIIIIIIIIIIM                                       
                       MIIIIIIIIIIIIIZZZZMMIIIIII7MMNMMIIIIZM                                      
                    MIIIIIIIIIIIIIIIINZZZOIIMIIIIIIIIIIIII7ZZIO                                    
                 .MIIIIIIIIIIIIIIIIIIMZZZZIIIINOIIIIIIIIIMZZZIIIM                                  

'
	
}

#####################################################################################################
