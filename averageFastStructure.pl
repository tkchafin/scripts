#!/usr/bin/perl 

use strict; 
use warnings;
use Getopt::Long;
use File::Basename; 
use Statistics::R; 

my @meanQ; 
my @log; 
my $out = "./avg_k";
my $help = 0; 
my $force =0; 
my $k; 
my $reps; 
parseArgs(); 

my ($filepath, $dirpath) = fileparse($meanQ[0]); 
@meanQ = glob("@meanQ"); 
@log = glob("@log"); 
my $count = 0; 
my @data; 
my $fnum = 0; 
foreach my $file(@meanQ){ 
  my @line; 
  if ($force == 0){
    $file !~ /.*meanQ/ and die "Error: File $file is missing .meanQ extension. Are you sure it is the correct file type? To skip this check, add the -f flag to your command-line call\n"; 
  }
  open (my $fh, $file) || die "Can't open $file\n";  
  $count++;
  my $lnum = 0; 
  while (<$fh>){ 
    $lnum++; 
    chomp; 
    @line = split /\s+/; 
    s/\s+//g; 
    next unless length; 
    if ($count==1){
      if (!defined $k){
        $k = @line unless $k;
        print "K value was not supplied; inferring clusters from file $file: $k\n"
      }
     } 
    @line != $k and die "Error: Line $lnum of file $file doesn't have the correct number of clusters ($k)\n";  
  } 
  #$fnum++; 
  close $fh; 
}

#Get likelihoods from .log files
my @lognames; 
my @likelihoods; 
foreach my $file (@log){
  my ($fpath, $dpath) = fileparse($file); 
  push @lognames, $dpath . $fpath; 
  open (my $fh, $file) || die "Can't open $file\n"; 
  while (<$fh>){
    chomp;
    if (m/Marginal Likelihood =/){ 
      s/Marginal Likelihood =//;
      push @likelihoods, $_;  
    }
  }
}

#Default use all reps if no subset number provided
if (!defined $reps){
  print "Warning: Number of replicates to subset not provided; using all by default\n";
  $reps = scalar(@lognames);
  #print $reps . "\n"; 
}

my $R = Statistics::R->new(); 
  $R->start; 
  $out = $out . $k . ".meanQ";
  $R->set('lognames', \@lognames); 
  $R->set('likelihoods', \@likelihoods); 
  $R->set('reps', $reps);
  $R->set('out', $out); 
  $R->send(q`options(scipen=999)`);
  $R->run(q`likes <- data.frame(lognames, likelihoods)`);
  $R->run(q`likes[,1] = sub(".log","",likes[,1])`); 
  #Set up R functions
  $R->send(q`
    ################################
    JSD.pair <- function(x, y){
	###Function to compute Shannon-Jensen Divergence
	###x and y are the frequencies for the same p categories
	u <- x/sum(x)
	v <- y/sum(y)
	m <- (u+v)/2
	if (all(u*v>0)){
		d <- (u*log(u/m)+v*log(v/m))/2
	} else {
		P1 <- u*log(u/m)
		P2 <- v*log(v/m)
		P1[is.nan(P1)] <- 0
		P2[is.nan(P2)] <- 0
		d <- (P1+P2)/2
	}
	return(sum(d))
    }
    ##############################
    matchPops=function(ga, gb, niter=3000) {
	### function to match population identifiers between fastStructure runs
	### based on permutations of column names and Shannon-Jensen divergences
	minsum=1000
	for (i in 1:niter) {
		names(gb)=sample(names(gb))
		sumjsd=0
		for (n in names(ga)) { 
			sumjsd=sumjsd+JSD.pair(ga[,n],gb[,n])
		}
		if (sumjsd<minsum) {
			minsum=sumjsd
			gbnames=names(gb)
		}
	}
	return(list("pops"=gbnames,"min.JSD"=minsum))
    }
    ##############################
    averageBest=function(likelihoods,top=25) {
	# matches populations assignments among best-likelihood runs,
	# averages assignemnt probabilities, returns averaged meanQ table
	bests=head(likes[order(likes[,2],decreasing=T),1],top)
	gs=read.table(paste(bests[1],".meanQ",sep=""))
	g1=gs
	print("top 1")
	for (b in 2:top) {
		gn=read.table(paste(bests[b],".meanQ",sep=""))
		names(gn)=matchPops(g1,gn)$pops
		gs=gs+gn[,names(g1)]
	}
	return(gs/top)
    }`
  );
  #Run averaging functions
  $R->run(q`means=averageBest(likelihoods=likes, top=reps)`);
  $R->run(q`write.table(means, file=out, sep="   ", quote=FALSE, na="NA", append=FALSE, row.names=FALSE,col.names=FALSE)`); 
 

#open (my $ofstream, ">$out") || die "Can't open $out\n"; 
#  for (my $i=0; $i<=$#data; $i++){ 
#    for (my $k=0; $k<=$#{$data[$i]}; $k++){ 
#      print $data[$i][$k]/$count . " ";
#    }
#    print "\n"; 
#  }


exit;

#########################################################################################

sub parseArgs{

my $message = 
"\n\nAverages multiple fastStructure runs for the same k value.  

If you have problems running the script let me know. It hasn't really been tested fully, and I threw it together quickly. 

Arguments

	-i	- Input fastStructure .meanQ files - wildcard usage is fine
	-o	- Output prefix and path
	-l	- Input fastStructure .log files
	-r	- Number of replicates to use. Script will choose top N reps based on likelihoods
	-k	- Provide a k value, otherwise it will be detected from column counts
	-f	- Shut up and stop checking files for .meanQ extension 
\n\n"; 

	my $result = GetOptions
	( 
	'i=s{1,}'	=> \@meanQ,
	'f!'		=> \$force,
	'l=s{1,}'	=> \@log, 
	'k=i'		=> \$k,
	'r=i'		=> \$reps,
	'o=s'		=> \$out, 
	'h!'		=> \$help
	);
@meanQ or die "\n\nNo meanQ specified!" . $message; 
@log or die "\n\nNo .log specified!" .  $message; 
$help == 1 and die $message; 

}


 
