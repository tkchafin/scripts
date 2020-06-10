#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#Initialize scalars
my $input;
my $popmap;
my $output="structure.in";
my $missing="-9";
my $suppress=0;
my $extra;
my $locnames=0;
my $popN = 1.0;
my $globalN = 1.0;
my $oneLine =0;
#Call sub parseArgs to parse command-line arguments
parseArgs();

#Some warnings...
if ($suppress == 0){
	$output eq "structure.in" and print "Warning: Output name not specified, using default of \"structure.in\"", "\n";
	$missing eq "-9" and print "Warning: Missing data value not given; using default of \"-9\"\n", "\n";
}

#Format output if default used
if ($output eq "structure.in"){
	my ($filepath, $dirpath) = fileparse ($input);
	$output = "$dirpath/$output";
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
my %popmap;
my %enum;
my %popcodes;
if ($popmap){
    open ( POPMAP, $popmap) || die "Derp: Can't open $popmap: $!";
    my $popcount = 1;
    while (<POPMAP>){
        chomp;
        my @c = split /\s+/, $_;
        if ($enum{$c[1]}){
            $popmap{$c[0]} = $enum{$c[1]};
						if (!exists $popcodes{$enum{$c[1]}}){
							$popcodes{$enum{$c[1]}} = $c[1];
						}
            #print "$c[0] is from pop# $popmap{$c[0]}\n";
        }else{
            $enum{$c[1]} = $popcount;
            $popcount++;
            $popmap{$c[0]} = $enum{$c[1]};
        }
    }
	close POPMAP;
	print "Population codes:\n";
		foreach my $p (sort keys %popcodes){
			print $p, ": ", $popcodes{$p}, "\n";
		}
}

#Begin going through phylip file
my $count = 0;

open ( OUTFILE, ">", $output) || die "Can't open $output: $!";
open ( PHY, $input ) || die "Can't open $input: $!";

my $samplecount = 0;
my $snpcount = 0;

#data structure to hold it, so we can print in order
my %structure;

while ( my $line = <PHY> ){
	$count++;
	$count == 1 and next; #Test if $count=1, if so then skip to next iteration

	#Split each line, store sequence name and sequence
	my @b = split /\s+/, $line;
	my @seq_array = split //, $b[1];

	#Build first line of structure file, containing "locus IDs"
	if ($count == 2){
	    my $locus_names= "\t\t";
	    for (my $i=1; $i <= scalar @seq_array; $i++){
		$locus_names .= "$i\t";
	    }
	    chop $locus_names;

	   if ($locnames == 1){
		 print OUTFILE "$locus_names\n";
	   }
	}
	#Begin building structure lines
	my $line_1 = "$b[0]\t";#Put in sequence name
	my $line_2 = "$b[0]\t";
	my $pop;

	if ($popmap){
		if (exists $popmap{$b[0]}){
			$pop = $popmap{$b[0]};
			#Add pop codes
			$line_1 .= "$popmap{$b[0]}\t";
			$line_2 .= "$popmap{$b[0]}\t";
		}
		else{
			next;
		}
	}

	if ($extra){
	   for (my $i=0; $i<$extra; $i++){
		$line_1 .= "0\t";
		$line_2 .= "0\t";
	    }
	}

	#Start adding allele data
	for( my $i=0; $i <= $#seq_array; $i++ ){
		if ($snpcount == 0){
			$snpcount = $#seq_array;
		}else{
			if ($snpcount != $#seq_array){
				print "Warning: Sample ",$b[0], " appears to have a different number of nucleotides. Something is wrong.\n";
			}
		}
		if ($first_line{ uc $seq_array[$i] }){
		    $line_1 .= "$first_line{ uc $seq_array[$i] }\t";
		}else{
		    $line_1 .= "-9\t";
		}
		if ($second_line{ uc $seq_array[$i] }){
		    if ($oneLine == 0){
		        $line_2 .= "$second_line{ uc $seq_array[$i] }\t";
		    }else{
			$line_1 .= "$second_line{ uc $seq_array[$i] }\t";
		    }
		}else{
		    if ($oneLine==0){
		   	$line_2 .= "-9\t";
		    }else{
			$line_1 .= "-9\t";
		    }
		}

	}

	chop $line_1;
	chop $line_2;

	if (not exists $structure{$pop}){
		$structure{$pop} = [];
	}

	if ($oneLine==0){
		#print OUTFILE $line_1, "\n";
		#print OUTFILE $line_2, "\n";
		push @{ $structure{$pop} }, [$line_1, $line_2];
	}else{
		push @{ $structure{$pop} }, [$line_1];
		#print OUTFILE $line_1, "\n";
	}
	$samplecount++;

}

foreach my $pop_key (sort {$a <=> $b} keys %structure) {
	foreach my $sample (@{$structure{$pop_key}}){
		foreach my $line (@{$sample}){
			print OUTFILE $line, "\n";
		}
	}
}

close PHY;
close OUTFILE;
print ("Done! Outputted ", $samplecount, " samples and ", $snpcount+1, " SNPs.\n");
exit;

############################SUBROUTINES######################################

sub parseArgs{
my $help=0;

my $usage= "\nUsage: $0 -i /path/to/phylip -p /path/to/popmap -o /path/to/output

The purpose of this script is to take a phylip-formatted file of concatenated SNPs (such as that output by the program pyRAD) and convert it to a structure-formatted file, with two lines for each individual representing the phased allele, as well as a column representing the a priori population/ locality assignment (as provided by the user in the form of a tab-delimited table).

Format of population map:
Sample1 	1
Sample2		1
Sample3		2
etc


Required Inputs

	-i, --input	-  Path to the input phylip file
	-p, --popmap	-  Path to the input population ID table
	-o, --output	-  Path to output (including desired filename)
	-n, --popN	- Percent missing data allowed per SNP per population [default=1.0]
			NOTE: Only applies when popmap provided
	-N, --globalN	- Percent missing data allowed per SNP globally [default=1.0]
			NOTE: N filters not implemented yet.

Optional inputs
	--oneLine	- Print phased alleles on one line
	-l, --loc	-  Bool, switch on printing of locus names in first row
	-e, --extra	-  Number of extra columns to insert
	-m, --missing	-  Desired code for missing data [Default is \"-9\"]
	-q, --quiet	-  Quiet mode; suppress internal warnings
	-x	- Exlude samples that are NOT in popmap

NOTE: Both gaps and N\'s will be coded as missing data.\n\n";

	my $result = GetOptions
	(
	'input|i=s'	=> \$input,
	'popmap|p=s'	=> \$popmap,
	'output|o=s'	=> \$output,
	'missing|m=s'	=> \$missing,
	'help|h!'	=> \$help,
	'extra|e=i'	=> \$extra,
	'loc|l!'	=> \$locnames,
	'quiet|q!'	=> \$suppress,
	'popN|n'	=> \$popN,
	'globalN|N' => \$globalN,
	'oneLine!' => \$oneLine
	);

	$help == 1 and die "$usage";
	$input || die "Input not specified!\n$usage";
};
