#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case no_auto_abbrev);
use Math::BigFloat;
use Data::Dumper;

my ($anno,$cuffdir,$cuffcompare,$tmap,$coverage,$cutoffstart,$cutoffend,$increment,$output,$chr) = ("","./","","trancripts.gtf.tmap","",0.01,0.20,0.01,"Optimizer","");
GetOptions (
  "a|anno:s" => \$anno,
  "cd|cuff_dir:s" => \$cuffdir,
  "cc|cuffcompare_prefix:s" => \$cuffcompare,
  "t|tmap_file:s" => \$tmap,
  "co|coverage_file:s" => \$coverage,
  "s|start:f" => \$cutoffstart,
  "e|end:f" => \$cutoffend,
  "i|increment:f" => \$increment,
  "o|output_directory:s" => \$output,
  "ch|chromosome:s" => \$chr,
);

$cutoffstart = Math::BigFloat->new($cutoffstart);###define floating point values at BigFloats to avoid embarrasing boundary issues.
$cutoffend = Math::BigFloat->new($cutoffend);
$increment = Math::BigFloat->new($increment);

##Lookup Cufflinks directory and check if cuffcompare prefix is valid, and the required files are present (transcripts.gtf, transcripts.gtf.tmap, CUFFCOMPARE.combined.gtf). Print error if any missing, this includes initial annotation file (-a) and the coverage file (-co).

###Check the Cufflinks directory and cuffcompare prefix are sensible.
if( $cuffdir =~ m/\/$/ ){
} else {
	$cuffdir = $cuffdir."/";
}

my $usage = "TSS-Optimizer.pl -a <Annotations.gtf> -cd <Cufflinks_directory> -cc <Cuffcompare_prefix> -t <tmap_file> -co <Coverage_file> -ch <chromosome>
Options:-
Short	Long			Default			Description
-s	start			[0.01]			X fraction of expected gene expression to begin trimming
-e	end			[0.20]			X fraction of expected gene expression to end trimming
-i	increment		[0.01]			Trimming increment
-o	output_directory	[./Optimizer]		Default output directory		
";

die "$usage\n" unless $anno && $cuffdir && $cuffcompare && $coverage && $chr;

if (-d $cuffdir ){
	if (-e  $cuffdir."transcripts.gtf"){
		if(-e $tmap ){
			if (! $cuffcompare eq ""){
				if(! -e $cuffdir.$cuffcompare.".combined.gtf"){
					die "Could not find $cuffcompare\.combined.gtf. Are you sure $cuffcompare is a valid cuffcompare prefix?";
				}
			} else {
				die "Please provide a cuffcompare prefix with -cc.";
			}
		}
		else{
			die "Could not find cuffcompare output $tmap, did you specify this file correctly?";
		}
	}
	else{
		die "Could not find either transcripts.gtf. Are you sure $cuffdir is a Cufflinks output directory???";
	}
} else {
	die "-cd Cufflinks directory does not exist";
} 

print "Cufflinks directory appears OK.\n";

###Check that coverage and annotation files exist.
if(! -e $anno){
	die "Please provide an annotation file with -a (This must be the one you used with Cuffcompare!)";
}

if(! -e $coverage){
	die "Please provide a coverage file (-co)";
}

##Create the output directory
$output = $output."-".$chr;
if( -d $output ){#Check the output directory exists, if it does then complain, otherwise make it.
	die "Directory $output already exists.";
} else {
	die "mkdir $output failed" unless (mkdir $output);
}

###Create a temporary coverage file with information for the chromosome of interest only
my $minicov = $output."\/".$chr."\.coverage";
print "Creating $minicov:- a coverage file for $chr\n";

print "$minicov\n";
open( TEMPCOV, ">$minicov" );
open( COV, "<$coverage" );
while( <COV> ){
	$_ =~ m/^(.*?)\t/;
	if( $1 eq $chr ){
		print TEMPCOV;
	}
}
close COV;
close TEMPCOV;

##Run Prep4Trim.pl (hard code the tmap_out and combined_out files)

my $tmap_out = "$output\/$chr.trim.gtf.tmap";
my $combined_out = "$output\/trim.combined.gtf";
my $prepcmd = "Prep4Trim.pl  -a ".$anno."  -i ".$cuffdir."transcripts.gtf  -t ".$tmap." -c ".$cuffdir.$cuffcompare.".combined.gtf  -to ".$tmap_out."  -go  ".$combined_out." -ch ".$chr;

print "Selecting Cufflinks predictions to trim.\n";
system( $prepcmd );


##Optimizer iterates over the main body of the code, incrementing $prefixstart by $iterator until $prefix > $cutoffend. Run Trim_Cufflinks => GetDistances for each value of $prefix.
while( $cutoffstart <= $cutoffend ){
##Run Trim_Cufflinks_GTF.pl 
	my $trimout = "$output\/$cutoffstart\.isotigs.trimmed.gtf";
	my $trimcmd = "Trim_Cufflinks_GTF.pl -t  ".$tmap_out." -g ".$combined_out." -ch ".$chr." -c ".$minicov." -f ".$cutoffstart." -o $trimout";
	print"Trimming Cufflinks annotations. Water level is $cutoffstart of expected per-base coverage for gene\n";
	system($trimcmd);

##Run Get-TSS.pl
	my $tssout = "$output\/$cutoffstart\.TSS.predictions.bed";
	my $tsscmd = "Get-TSS.pl -t ".$tmap_out." -g ".$trimout." -o ".$tssout;
	print "Getting transcription start sites for $cutoffstart cutoff.\n";
	system( $tsscmd );

##Run GetDistances.pl 
	my $distout = "$output\/$cutoffstart\.distances.txt";
	my $distcmd = "GetDistances.pl -a ".$anno." -p ".$tssout." -o ".$distout;
	print "Calculating distances between predicted and nearest annotated TSS for trimming at $cutoffstart.\n";
	system( $distcmd );

	$cutoffstart = $cutoffstart + $increment;
}


##CollateDistances.pl collates the distances files and produces a melted data frame (ggplot2 style) for all optimizer iterations.

my $collatecmd = "CollateDists.pl $output $output\/*.distances.txt";
system ( $collatecmd );


##Pass the CutoffDistances.txt file made by CollateDistances.pl into R script to get MAD plot
chdir $output;
system( "MADPlot.R" );

exit;
