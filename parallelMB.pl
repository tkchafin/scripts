#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Cwd;
use File::Basename; 
use Parallel::ForkManager; 

#runMRBAYES.pl runs mrbayes for a given gene. It will also write a mrbayes command block to the nexus file, if it isn't already present, and run mbsum to combine independent runs from a gene as well as to calculate P(T|gene) 
#Requires that mrbayes and mbsum be in your PATH.

# Define defaults for variables
	my @input;
	my $procs=1; 
	my $batch="1";
	my $ngen="1000000";
	my $nruns="2";
	my $nst="6";
	my $rates="Invgamma";
	my $savebrlens="yes";
	my $printfreq="1000";
	my $samplefreq="100";
	my $nchains="4";
	my $burnin;
	my $relburnin="yes";
	my $burninfrac=".25";
	my $cmd = "mb"; 
	my $BURNIN; 
	my $write_block=0; #If 0 then write NEXUS block
	my $nexus=0; #If 1 then convert files from FASTA to NEXUS
	my $workdir; 
	my $outdir; 
	my $help=0; 
	my $no_mbsum=0; 
	my $cleanup=0;
	my $sumburn=0; 
	
parseArgs();	#Subroutine to call GetOptions

#	print "infile = $input\n";
#	print "batch = $batch\n";
#	print "ngen = $ngen\n";
#	print "savebrlens = $savebrlens\n";
#	print "printfreq = $printfreq\n";
#	print "samplefreq = $samplefreq\n";
#	print "nchains = $nchains\n";
#	print "relburnin = $relburnin\n";
#	print "rates = $rates\n";
#	print "nst = $nst\n";
	
#Get file base name and dirpath
my ($filename, $dirpath) = fileparse ($input[0]);
 
#Format vars
$ngen = "ngen=$ngen";
$nst = "nst=$nst";
$rates = "rates=$rates";
$savebrlens = "savebrlens=$savebrlens";
$printfreq = "printfreq=$printfreq";
$samplefreq = "samplefreq=$samplefreq";
$nchains = "nchains=$nchains";
$relburnin =~ /^n/i and $BURNIN = "relburnin=no burnin=$burnin";
$relburnin =~ /^y/i and $BURNIN = "burninfrac=$burninfrac";
$nruns = "nruns=$nruns"; 
$workdir or $outdir = "$dirpath";
$workdir and $outdir = "$workdir";
$sumburn = "-n $sumburn";  
@input = glob "@input"; 

#Capture file number/ name
 
mkdir "$outdir/mb.$batch" unless ( -d "$outdir/mb$batch" );
print "\n";  

#Check if filetype = Nexus, if not then convert


#######################MrBayes Independently on each locus#######################

#Format some variables
if ($cleanup == 0){ $cleanup = "";
	}else{ $cleanup = "-k";
}

         mkdir "$outdir/mb.$batch";

#Format to NEXUS if necessary

$nexus ==1 and print "\nparallelMB.pl: Converting files to NEXUS format...\n";
$nexus == 1 and system("fasta2nexus.pl -i @input -p $procs") and print "Done!\n\n";

#print "Running mrbayes on each locus.........\n";
my $ID; 

my $pm = Parallel::ForkManager->new($procs);
	
     my @nex = glob "$dirpath/*.nex";

     foreach my $file (@nex){
	
	my $pid = $pm-> start and next; 

	my ($filename, $dirpath) = fileparse ($file);
	MB_block ( $file ) unless $write_block =~ "1";

	$filename =~ /(\w+).\w/;
	my $ID = "$1";
	
	print "parallelMB.pl: Submitting $file to MRBAYES...\n";
        system("$cmd $file > $outdir/mb.$batch/mb.$ID.batch$batch.out");

     $pm->finish;
     }

$pm->wait_all_children;
    

#######################Post processing###########################
#Some clean up of files we don't need if "cleanup" option was toggled 

chdir ( $dirpath ); 

if ($cleanup == 1){
system("rm -f $filename.*.p"); 
system("rm -f $filename.ckp*");
unlink("$filename.mcmc");          
}

`mv -f $filename.* $outdir/mb.$batch/.` ;  #Move files to keep to output dir



#########################mbsum####################################
 
if ($no_mbsum != "1"){  #Check if no_mbsum was toggled on

chdir "$outdir/mb.$batch"; #Change directories for less typing in next few lines...

#mbsum system call
system ( "mbsum $sumburn -o $filename.in $filename*.t >> mb.$ID.batch$batch.out");
$cleanup == 1 and system ("rm -f $outdir/mb.$batch/$filename*.t");
}






exit; 

########################################################################
###########################SUBROUTINES###################################
########################################################################

sub parseArgs{ 
	#Message to print if mandatory variables not declared
	my $usage ="\nUsage: $0 --input=/path/to/*.nex --batch=batch_number --option=value [options]  

Type --help, or -h to view all options. 


Mandatory 
	-i, --input		-  path to the input files in nexus format 
	-p, --procs		-  Number of processes for parallel execution	

General Options
	--batch			-  batch number must be specified, if different than default [Default=1]
	-c, --cmd		-  command to call mrbayes, if different than default [default=mb]
	-o, --out		-  Specify an output directory to store the MrBayes outputs [default=same folder as nexus file]
	-x, --skip		- Indicates that MrBayes command block is already present in NEXUS files. If this option is toggled on, then a new block will not be written [Usage: --skip, or -x] 
	-e, --nex		- Indicates that files need to be converted to NEXUS format (from FASTA) Usage: --nex or -e		- 
	-k, --cleanup		- Tells runMRBAYES.pl to cleanup files which are not needed for BUCKy Bayesian Concordance Analyses (.p files, .mcmc, .ckp, and .t files after summing) [Usage: --cleanup, or -k] 
	--no_mbsum		- If toggled on, will skip mbsum step and just output .t files
	
lset options - Sets likelihood model parameters - Change if different than default
	-m, --nst 		- sets the substitution model [must be equal to 2, 4, or 6; Default=6(GTR)]
	-t, --rates		- options for substitution rates: equal, gamma, propinv, invgamma, adgamma [Default=Invgamma] 

mcmcp options - Sets MCMC paramenters - Change if different than default
	-n, --ngen 		- number of generations for each MCMC chain to run [default=1000000]
	-r, --nruns		- number of independent analyses to run [default=2] 
	--savebrlens 		- do you want to save the branch lengths?  Options: yes or no
	--printfreq 		- how often do you want to print the likelihood values? [must be integer; default=100]
	-s, --samplefreq 	- how often should the likelihood values be samples? [must be integer; default=100]
	-c, --nchains		- how many independent MCMC chains do you want to run? [default=4]
	--relburnin 		- do you want to discard burnin values as a fraction?  options: yes or no [default=yes]
	-b, --burninfrac 	- if relburnin=yes what fraction of the MCMC chain should be discarded as burnin [default=0.25]
	--burnin		- if relburnin=no how many samples to discard as burnin

mbsum options - Combines all runs from a gene and calculates P(T|gene) with a standardized representation of each tree 

	-u, --sumburnin		- Number of trees for mbsum to discard (burn-in) [default=0] 

\n";

	my $options = GetOptions 
		( 
		'input|i=s{1,}'		=>	\@input,
		'batch=s'		=>	\$batch,
		'out|o=s'		=>	\$workdir,
		'cmd|c=s'		=>  	\$cmd,
		'ngen|n=s'		=>	\$ngen, # how many MCMC generations to run
		'nst|m=s'		=>	\$nst, # set the substitution model
		'rates|t=s'		=>	\$rates, #set the model for substitution rates
		'nruns|r=s'		=> 	\$nruns,
		'savebrlens=s'		=>	\$savebrlens, # Do you want to save the branch lengths?
		'printfreq=s'		=>	\$printfreq, # How often to print to the file
		'samplefreq=s'		=>	\$samplefreq, # How often should likelihood values be sampled
		'nchains|c=s'		=>	\$nchains, # How many independent chains to run
		'burnin=s'		=>	\$burnin, # how many samples should be discarded as burnin
		'relburnin=s'		=>	\$relburnin, # should burnin samples be discarded as fraction?  Options: yes or no
		'help|h!'		=>	\$help,
		'no_mbsum!'		=>	\$no_mbsum,
		'cleanup|k!'		=>	\$cleanup,
		'sumburnin|u=i'		=>	\$sumburn,
		'skip|x!'		=> 	\$write_block,
		'nex|e!'		=> 	\$nexus,
		'procs|p=i'		=> 	\$procs, 
		'burninfrac|b=s'	=>	\$burninfrac, # what fraction of the MCMC chain should be discarded as burnin
		);
	
	$help =~ "1" and die "\n$usage\n";
	@input or die "\n\nDerp: Input not specified!\n\n$usage\n"; 
	if ( $input[0] !~ /.nex$/){ $nexus == 0 and die "/n/nDerp: Infiles do not have the extension .nex, are you sure you don't need to toggle \"--nex\"?\n$usage\n"}; 
	
}
#################################################################################################

#Writes a mrbayes command block to nexus file
sub MB_block { 

my $file = $_[0];

print "parallelMB.pl: Adding command block to $file.\n";

my ($filename, $dirpath) = fileparse ($file);

open ( NEX, ">>$file" ) || die "Derp: Can't open $file for writing!\n";

print NEX "begin mrbayes;       
        set autoclose=yes nowarn=yes quitonerror=no;
      	lset $nst $rates;
        mcmcp $ngen $printfreq $nruns $BURNIN $savebrlens $nruns $samplefreq $nchains;
	mcmc;
	sump; 
	sumt; 
end;\n";

};
