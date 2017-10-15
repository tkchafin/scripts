#!/usr/bin/perl 

use strict; 
use warnings; 
use Getopt::Long qw( :config posix_default no_ignore_case );
use File::Basename;

#Initialize scalars
my $input; 
my $popmap;
my $output="structure.in";
my $missing="-9";  
my $suppress=0; 
my $type; 
my $snp = 0; 
#Call sub parseArgs to parse command-line arguments
parseArgs();


my $message = "#File created by seq2structure.pl; script by Tyler K. Chafin last updated 12-Dec-14";

#Some warnings...
if ($suppress == 0){
	$output eq "structure.in" and print "Warning: Output name not specified, using default of \"structure.in\"", "\n";
	$missing eq "-9" and print "Warning: Missing data value not given; using default of \"-9\"\n", "\n";
}

#Format output if default used
if ($output eq "structure.in"){
	my ($filepath, $dirpath) = fileparse ($input); 
	#$output = "$output"; 
};

#Specify iupac abiguity codes and how to write them out to structure file
my $iupac="A       1 1
C       2 2
G       3 3
T       4 4
N       $missing $missing
-       $missing $missing
R       1 3
Y       2 4
S       2 3
W       1 4
K       3 4
M       1 2 ";


my %first_line; 
my %second_line;

#Build hashes of above iupac codes
for my $line (split "\n", $iupac){
	chomp $line; 
	my @a = split /\s+/, $line; 
	$first_line{ $a[0] } = $a[1]; 
	$second_line{ $a[0] } = $a[2];

};

#Store population identifiers for each individual (from popmap)
open ( POPMAP, $popmap) || die "Derp: Can't open $popmap: $!";
	my %popmap;
	    while (<POPMAP>){
	    chomp;
	    my @c = split /\s+/, $_;
	    $popmap{$c[0]} = $c[1];
	#print "$c[0] is from pop# $popmap{$c[0]}\n";
}
close POPMAP;

#Begin going through phylip file
my $count = 0;
my @b; 
my @seq_array;

open ( OUTFILE, ">$output") || die "Can't open $output: $!";
open ( PHY, $input ) || die "Can't open $input: $!";

while ( my $line = <PHY> ){
	$count++; 
		
	chomp $line; 
	#Split each line, store sequence name and sequence 
	if ($type =~ /p/i){
		$count == 1 and next; #Test if $count=1, if so then skip to next iteration.
		@b = split /\s+/, $line; 
		@seq_array = split //, $b[1]; 
	
		#Build first line of structure file, containing "locus IDs"
		if ($count == 2){
			#my $locus_names= "\t\t";
			#for (my $i=1; $i <= scalar @seq_array; $i++){
				#$locus_names .= "$i\t"; 
			#}
			#chop $locus_names;
			#print OUTFILE "$message\n";
			#print OUTFILE "$locus_names\n";
		}
	}
	
	if ($type =~ /f/i){
		
		if ($count ==2){
			#my $locus_names= "\t\t";
		
			#for (my $i=1; $i <= scalar(@seq_array); $i++){
				#$locus_names .= "$i\t"; 
			#}
			#chop $locus_names;
			#print OUTFILE "$message\n";
			#print OUTFILE "$locus_names\n";
		}
		
		if ($line =~ /^\>(\S+)/){ 
			$b[0] = $1; 
			next; 
		}elsif ($line =~ /[ACGT]+/i){ 
			@seq_array = split //, $line; 
		}else{
			
			next; 
		}
		
	}
		
		 
	#Begin building structure lines
	my $line_1 = "$b[0]\t";#Put in sequence name
	my $line_2 = "$b[0]\t";

	#Add pop codes
	$line_1 .= "$popmap{$b[0]}\t";
	$line_2 .= "$popmap{$b[0]}\t";

	#Start adding allele data
	for( my $i=0; $i <= $#seq_array; $i++ ){ 
		$line_1 .= "$first_line{ uc $seq_array[$i] }\t"; 
		$line_2 .= "$second_line{ uc $seq_array[$i] }\t";  
	}
	
	chop $line_1; 
	chop $line_2;
	
	print OUTFILE $line_1, "\n"; 
	print OUTFILE $line_2, "\n"; 
	print "Sample $b[0] done...\n";
	$count++;
}

close PHY;
close OUTFILE;

#If SNP check toggled on, rewrite file with only snps
my $loci=0;
if ($snp == 1){
	open (STR, "$output") || die "Cannot open $output for reading: $!\n";
	my $comments = "";
	my @data;  
	my $num=0; 
	my $locnames;
	foreach (<STR>){ 
		chomp; 
		#If line is a comment, capture to reprint later
		$_ =~ /^#/ and $comments .= $_ and next; 
		#If line contains variable number of spaces and nothing else, skip
		$_ =~ /^ *$/ and next; 
		#If column has locus names
		#$num == 0 and $locnames = $_; 
		#Capture elements in line
		my @line = split("\t"); 
		
		#Build array of arrays with secondary arrays as the columns from structure file
		for (my $col = 0; $col < scalar(@line); $col++){ 
			 push(@{$data[$col]}, $line[$col]); 
		}
		$num++; 
	}
	close STR;
	 
	# Check each column for unique 	
	for (my $col = 0; $col < scalar(@data); $col++){ 
		$col < 2 and next; #Skip sample and popID columns
		my %counts; 
		$counts{$_}++ for @{$data[$col]};
			#print keys(%counts) ."\n";
		my $number = keys %counts; 
		#If column doesn't contain a SNP, delete it.
		unless ($number > 1){
			undef $data[$col];
			next;
		}
		$loci++;
	}
	
	#Build new structure file containing only the SNPs
	my $ind = (scalar(@{$data[0]})/2);
	print "\n######################################\n\n";
	print "Number of Individuals: $ind\nNumber of SNPs discovered: $loci\n";
	my $outfile = "N" . $ind ."-" . "L" . $loci . "_" . "$output";
	print "\nWriting $outfile...\n\n";
	open (NEWOUT, ">$outfile") || die "Can't open $output for re-writing: $!\n";
	$comments and print NEWOUT "$comments\n"; 
	#print NEWOUT "$locnames\n";
	for (my $row = 0; $row < scalar(@{$data[1]}); $row++){ 
		for (my $col = 0; $col < scalar(@data); $col++){ 			
			if (defined $data[$col][$row]){
				print NEWOUT $data[$col][$row] . "\t";
			}
		}
		print NEWOUT "\n"; 
	}
	close NEWOUT;
}
	

exit;  

############################SUBROUTINES######################################

sub parseArgs{

my $help=0;

my $usage= "\nUsage: $0 -i /path/to/seqfile -p /path/to/popmap -o /path/to/output

The purpose of this script is to take a phylip or fasta-formatted file of concatenated SNPs (such as that output by the program pyRAD) and convert it to a structure-formatted file, with two lines for each individual representing the phased allele, as well as a column representing the a priori population/ locality assignment (as provided by the user in the form of a tab-delimited table).

Format of population map: 
Sample1 	1
Sample2		1
Sample3		2
etc


Required Inputs

	-i, --input	-  Path to the input sequencefile
	-p, --popmap	-  Path to the input population ID table
	-o, --output	-  Path to output (including desired filename)
	-t, --type	-  Input file type (phylip or fasta)

Optional inputs

	-m, --missing	-  Desired code for missing data [Default is \"-9\"]
	-s, --snp	-  Check for SNPs; only write snps to str file
	-q, --quiet	-  Quiet mode; suppress internal warnings

NOTE: Both gaps and N\'s will be coded as missing data.
NOTE: Script assumes a perfect alignment (same length, gaps and N's inserted where needed).
NOTE: SNP checking currently not functional. 
NOTE: Script does not create a row in structure file for locus names. Will add this functionality back in later, if necessary.
TODO: Add built-in check for filetype so it doesn\'t need to be specified.\n\n";

	my $result = GetOptions
	(
	'input|i=s'	=> \$input,
	'popmap|p=s'	=> \$popmap,
	'output|o=s'	=> \$output,
	'missing|m=s'	=> \$missing,
	'help|h!'	=> \$help,
	'snp|s!'	=> \$snp, 
	'type|t=s'	=> \$type,  
	'quiet|q!'	=> \$suppress,
	);

	$help == 1 and die "$usage";
	$input || die "Input not specified!\n$usage";
	$popmap || die "Popmap not provided!\n$usage";
	$type || die "Popmap not provided!\n$usage";
};
	
