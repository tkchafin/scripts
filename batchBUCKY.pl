#!/usr/bin/perl

use strict; 
use warnings; 
use Getopt::Long; 
use File::Basename; 

our $cmd="bucky"; 
our $input; 
our $ngen=100000; 
our $nrun=2;
our $alpha=1;  
our $help=0; 
our $nchain=1; 
our $rate="";
our $alpham=""; 
our @other = ();
our $cutoff="0.05";
our $ind=0;
our $spacesaver=0;
our $out="bca";

#Parse command line arguments to define above variables

parseArgs(); 

#Format some variables

my $other = "@other";
$rate =~ /\S/ and $rate = "-r $rate";
$alpham =~ /\S/ and $alpham = "-m $alpham"; 
$ind == 1 and $ind = "--use-independence-prior";
$ind == 0 and $ind = ""; 
$spacesaver == 1 and $spacesaver = "--opt-space";
$spacesaver == 0 and $spacesaver = "";
my ($filepath, $dirpath) = fileparse($input);
$input = "$dirpath\*.in";
my ($outname, $outpath) = fileparse($out); 

#BUCKy system call 
print $outpath, "\n"; 
chdir "$outpath"; 

system ("$cmd -a $alpha -k $nrun -c $nchain $rate $alpham -o $outname -cf $cutoff $ind $other $spacesaver $input");




exit;

#####################################SUBROUTINES##############################################

sub parseArgs{

    my $usage="\nUsage: $0 -i /path/to/*.in [-option value] or [--option=value] 

batchBUCKY.pl takes as input the summed .t files from mrbayes (summarized per locus via mbsum) and performs a Bayesian Concordance Analysis to assess what proportion of the genome supports different phylogenetic topologies. 

--------------------------------------Mandatory Input-------------------------------
	-i, --input	- Path to .in files created by mbsum (automatically generated by runMRBAYES.pl)


---------------------------------------General Options-------------------------------
	--cmd		- Command to call bucky, if different than default [default=bucky] 
	-o, --out		- Output file root name [Default=bca] Can also include path to output directory [e.g. -o /path/to/bca]
	
----------------------------------------BUCKy Options-------------------------------- 
	-a, --alpha	- Use this option to set the a priori level of discordance among loci [default=1] 
	-n, --ngen	- Number of generations for MCMC. Burnin will automatically be 10% of the desired number of post-burnin updates. [default=100,000] 
	-k, --nrun	- Number of independent analyses to run 
	-f, --cutoff	- Provide a cutoff Concordance Factor value, above which all splits will be retained [default=0.05] 
	--other		- Use this option to set any other bucky parameters or functions [example: --other -s1 1234 --calculate-pairs --create-single-file] 	
	--ind		- Use independent priors. Assumes a priori that loci have independent histories. [Usage: --ind ] 
	--opt-space	- This option accomodates large data sets with space optimization. [Usage: --opt-space ] 

----------------------------------------MCMCMC Options--------------------------------
	-c, --nchain	- This option toggles on Metropolic coupled MCMC. Any chains more than one will be \"hot\" chains which will occasionally swap states with the cold chain to improve mixing [Default=1; i.e. no heated chains]
	-r, --rate	- If MCMCMC is used, this controls the rate at which chains swap [default=100]
	-m, --alpham	- Heated chains in MCMCMC use higher alpha values than the cold chain. This parameter sets the multiplier for the heated alpha value [default=10]\n\n"; 

	my $result = GetOptions 
	(	
	'input|i=s'	=> \$input,
	'cmd=s'		=> \$cmd, 
	'alpha|a=s'	=> \$alpha,
	'ngen|n=i'	=> \$ngen, 
	'nrun|k=i'	=> \$nrun,
	'cutoff|f=s'	=> \$cutoff,
	'other=s{1,}'	=> \@other,
	'nchain|c=i'	=> \$nchain,
	'rate|r=i'	=> \$rate,
	'alpham|m=i'	=> \$alpham,
	'help|h!'	=> \$help,
	'out|o=s'	=> \$out,
	'ind!'		=> \$ind,
	'opt-space!'	=> \$spacesaver,
	); 

$help == 1 and die "\n$usage\n";
$input or die "\nInput not specified!\n$usage\n";

}
