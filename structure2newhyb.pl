#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

my %opts;
getopts('s:o:a:b:x:r:h', \%opts);

if ($opts{h}){
	&help;
	die "Exiting because help menu was called.\n\n"
}

my ($str, $out, $A, $B, $X, $rand) = &parseArgs(\%opts);

#Arrays of pop IDs for each class
my @popA = split /\+/, $A;
my %mapA = map {$_ => 1} @popA;
my @popB = split /\+/, $B;
my %mapB = map {$_ => 1} @popB;
my @popX = split /\+/, $X;
my %mapX = map {$_ => 1} @popX;

print "Parental population 1: ";
foreach (@popA){
	print $_, " ";
}
print "\n";

print "Parental population 2: ";
foreach (@popB){
	print $_, " ";
}
print "\n";

print "Admixed population(s): ";
foreach (@popX){
	print $_, " ";
}
print "\n\n";

#Instantiate hashes to hold the data
my %dataA;
my %dataB;
my %dataX;
my $numLoci;

print "Parsing structure file... ";

open (my $fhs, $str) or die "Can't open $str\n";

my @line1 = ();
my @line2 = ();
while (my $row = <$fhs>){
	chomp $row;
	if (!@line1){
		@line1 = split "\t", $row;
		if (!$numLoci){
			$numLoci = (scalar(@line1)-2);
		}else{
			if ($numLoci != (scalar(@line1)-2)){
				die "Number of columns in each row not the same.\n"
			}
		}
		next;
	}elsif(!@line2){
		@line2 = split "\t", $row;
		if ($numLoci != (scalar(@line2)-2)){
			die "Number of columns in each row not the same.\n"
		}
		if ($line1[0] ne $line2[0]){
			die "$line1[0] not equal to $line2[0] something is wrong!\n";
		}
	}else{
		@line1 = ();
		@line2 = ();
		@line1 = split "\t", $row;
		next;
	}

	my $sample = $line2[0];
	my $pop = $line2[1];

	if (exists($mapA{$pop})){
		$dataA{$sample} = &fetchNewHybLine(\@line1, \@line2);
	}elsif (exists($mapB{$pop})){
		$dataB{$sample} = &fetchNewHybLine(\@line1, \@line2);
	}elsif (exists($mapX{$pop})){
		$dataX{$sample} = &fetchNewHybLine(\@line1, \@line2);
	}else{
		@line1 = ();
		@line2 = ();
		next;
	}
}
close $fhs;
print ("Done.\n");

#Check that all hashes have data.
%dataA or die "Error: No data for parental population 1!\n";
%dataB or die "Error: No data for parental population 2!\n";
%dataX or die "Error: No data for admixed population!\n";

print ("Removing loci with no data for any population... ");
my @blacklist;
for (my $k = 0; $k < $numLoci; $k++){
	#Population A
	my $bad = 1;
	for my $key (keys %dataA){
		if ($dataA{$key}->[$k] != 0){
			$bad = 0;
		}
	}
	if ($bad == 1){
		push @blacklist, $k;
	}

	#Population B
	$bad = 1;
	for my $key (keys %dataB){
		if ($dataB{$key}->[$k] != 0){
			$bad = 0;
		}
	}
	if ($bad == 1){
		push @blacklist, $k;
	}

	#Population X
	$bad = 1;
	for my $key (keys %dataX){
		if ($dataX{$key}->[$k] != 0){
			$bad = 0;
		}
	}
	if ($bad == 1){
		push @blacklist, $k;
	}
}

#Removing monomorphic loci
print "Removing monomorphic loci... ";
my $countM = 0;
for (my $k = 0; $k < $numLoci; $k++){
	#Population A
	my $var = 0;
	my $obs = 0;
	for my $key (keys %dataA){
		if ($obs != 0){
			$obs = $dataA{$key}->[$k];
		}
		if ($obs != $dataA{$key}->[$k]){
			$var = 1;
		}
	}
	if ($var == 1){
		next;
	}

	#Population B
	for my $key (keys %dataB){
		if ($obs != 0){
			$obs = $dataB{$key}->[$k];
		}
		if ($obs != $dataB{$key}->[$k]){
			$var = 1;
		}
	}
	if ($var == 1){
		next;
	}

	#Population X
	for my $key (keys %dataX){
		if ($obs != 0){
			$obs = $dataX{$key}->[$k];
		}
		if ($obs != $dataX{$key}->[$k]){
			$var = 1;
		}
	}
	if ($var == 0){
		$countM++;
		push @blacklist, $k;
	}
}
print "Done.\n";

my %seen;
@seen{@blacklist} = ();
my $numLost = scalar(keys %seen);
print "$numLost loci total loci dropped for missing data, or lacking variation.\n";

#Check if random sampling selected
my %skip;
if ($rand){
	my $tot = $numLoci - $numLost;
	if ($tot < $rand){
		print "Warning: You chose to randomly sample $rand loci, but only $tot loci remain. Skipping random sampling.\n";
		%skip = %seen;
	}else{
		print "Randomly sampling $rand loci... ";
		%skip = randomSampleLoci($numLoci, \%seen, $rand);
		my $len = scalar(keys %skip);
		#print "There are $len loci to skip of $numLoci total...\n";
		print "Done.\n";
	}
}else{
	%skip = %seen;
}

print ("Writing NewHybrids outputs with prefix \"$out\"... ");

my $ind = $out . "_individuals.txt";
my $dat = $out . ".txt";
open (my $ifh, ">$ind") or die "Can't open $ind for writing\n";
open (my $dfh, ">$dat") or die "Can't open $dat for writing\n";

my $numInds = scalar(keys %dataA) + scalar(keys %dataB) + scalar (keys %dataX);

#Print header to data file
my $keptLoci = $numLoci - scalar(keys %skip);
print $dfh "NumIndivs $numInds\n";
print $dfh "NumLoci $keptLoci\n";
print $dfh "Digits 1\n";
print $dfh "Format Lumped\n";
print $dfh "LocusNames";
for (my $i = 0; $i < $numLoci; $i++){
	print $dfh " Locus$i" unless exists $skip{$i};
}
print $dfh "\n";

#Print values for Pop A
my $indCount = 1;
for my $key (keys %dataA){
	print $ifh "Pure_PopA_$key\n";
	print $dfh "$indCount";
	my $index = 0;
	for my $loc (@{$dataA{$key}}){
		print $dfh " $loc" unless exists $skip{$index};
		$index++;
	}
	print $dfh "\n";
	$indCount++;
}

#Values for pop B
for my $key (keys %dataB){
	print $ifh "Pure_PopB_$key\n";
	print $dfh "$indCount";
	my $index = 0;
	for my $loc (@{$dataB{$key}}){
		print $dfh " $loc" unless exists $skip{$index};
		$index++;
	}
	print $dfh "\n";
	$indCount++;
}

#Values for pop X
for my $key (keys %dataX){
	print $ifh "Admix_$key\n";
	print $dfh "$indCount";
	my $index = 0;
	for my $loc (@{$dataX{$key}}){
		print $dfh " $loc" unless exists $skip{$index};
		$index++;
	}
	print $dfh "\n";
	$indCount++;
}

close $ifh;
close $dfh;

print "Done.\n\n";
exit;

###############################################################################
################################ Subroutines ##################################
###############################################################################

# subroutine to print help
sub help{

  print "\nPerl script to convert 2-line format Structure file to inputs for NewHybrids. Structure format should have no extra columns or rows, and 1 column with a population identifier.\n\n";
	print "Author: Tyler K. Chafin\n";
	print "Program Options:\n";

  print "\t-s:\tTab-delimited Structure file\n";
  print "\t-o:\tPrefix for output files [default=newhyb]\n";
	print "\t-a:\tInteger identifier(s) for Parent population 1 [Multiple as: -a 1+2+3]\n";
	print "\t-b:\tInteger identifier(s) for Parent population 2 [Multiple as: -b 1+2+3]\n";
	print "\t-x:\tInteger identifier(s) for admixed population [Multiple as: -x 1+2+3]\n";
	print "\t-r:\tInteger to randomly sample from loci [default=not used]\n";
	print "\t-h:\tBoolean. Calls help menu.\n\n";

}

# subroutine to parse the command line options
sub parseArgs{

  my( $params ) =  @_;
  my %opts = %$params;

	my $s = $opts{s} or die "Structure file not given\n";
	my $o = $opts{o} || "newhyb";
	my $a = $opts{a} or die "Identifiers for parent 1 not given\n";
	my $b = $opts{b} or die "Identifiers for parent 2 not given\n";
	my $r = $opts{r} || "";
	my $x = $opts{x} or die "Identifiers for admixed pop not given\n";

  return( $s, $o, $a, $b, $x, $r);

}

# subroutine to parse the command line options
sub fetchNewHybLine{
	my $line1 = $_[0];
	my $line2 = $_[1];
	my @export;
	my $counter = 0;
	for my $col (@{$line1}){
		if ($counter < 2){
			$counter++;
			next;
		}
		my $col_s = NHsanitize($col);
		my $other = $line2->[$counter];
		my $other_s = NHsanitize($other);
		#$other =~ s/-9/0/g;
		my $nh = $col_s . $other_s;
		push @export, $nh;
		$counter++;
	}
	return (\@export);
}

sub NHsanitize{
	my $s = $_[0];
	my $new = $s;
	$new =~ s/0/\!/g;
	$new =~ s/1/\@/g;
	$new =~ s/2/\#/g;
	$new =~ s/3/\$/g;
	$new =~ s/-9/0/g;
	$new =~ s/\!/1/g;
	$new =~ s/\@/2/g;
	$new =~ s/\#/3/g;
	$new =~ s/\$/4/g;
	return ($new);

}


# subroutine to randomly sample X loci to keep, and return new blacklist
sub randomSampleLoci{
	my $total = $_[0];
	my $blacklist = $_[1];
	my $sample = $_[2];
	my %new_blacklist = %$blacklist;

	my %kept;
	my $pick;

	while (scalar(keys %kept) < $sample){
		$pick = int(rand($total))+1;
		if (!$new_blacklist{$pick}){
			$kept{$pick} = undef;
		}
	}
	for (my $i=1;$i<=$total; $i++){
		if (!exists $kept{$i}){
			$new_blacklist{$i} = undef;
		}
	}
	return (%new_blacklist);
}
