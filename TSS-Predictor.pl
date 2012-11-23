#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case no_auto_abbrev);
use Data::Dumper;

my ($anno,$cuffdir,$cuffcompare,$tmap,$coverage,$cutoff,$upstream,$downstream) = ("","./","","trancripts.gtf.tmap","",0.1,0,1);
GetOptions (
  "a|anno:s" => \$anno,
  "cd|cuff_dir:s" => \$cuffdir,
  "cc|cuffcompare_prefix:s" => \$cuffcompare,
  "t|tmap_file:s" => \$tmap,
  "co|coverage_file:s" => \$coverage,
  "f|cutoff_fraction:s" => \$cutoff,
  "u|upstream_bases:s" => \$upstream,
  "d|downstream_bases:s" => \$downstream,
);

my $tmap_out = "trim.gtf.tmap";
my $combined_out = "trim.combined.gtf";
my $predictions_out = "tss.prediction.bed";
my $distances_out = "distances.txt";


##Lookup Cufflinks directory and check if cuffcompare prefix is valid, and the required files are present (transcripts.gtf, transcripts.gtf.tmap, CUFFCOMPARE.combined.gtf). Print error if any missing, this includes initial annotation file (-a) and the coverrage file (-co).

###Check the Cufflinks directory and cuffcompare prefix are sensible. Document my code wth POD???
if( $cuffdir =~ m/\/$/ ){
} else {
	$cuffdir = $cuffdir."/";
}

my $usage = "TSS-Predictor.pl -a <Annotations.gtf> -cd <Cufflinks_directory> -cc <Cuffcompare_prefix> -t <tmap_file> -co <Coverage_file>
Options:-
Short	Long			Default			Description
-f	cutoff_fraction		[0.1]			Fraction of expected coverage to trim annotation
-u	upstream_bases		[1]			Number of bases upstream of predicted TSS to report
-d	downstream_bases	[0]			Number of bases downstream of predicted TSS to report
-po	predictions_out		[tss.prediction.bed]	TSS predictions file
-do	distances_out		[distances.txt]		Distances from predicted transcription start sites to nearest annotated transcription start site
";

die "$usage\n" unless $anno && $cuffdir && $cuffcompare && $coverage;

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

##Get current working directory
#my $wd = `pwd`;
#$wd =~ s/\n/\//;

##Run Prep4Trim.pl
my $prepcmd = "Prep4Trim.pl  -a ".$anno."  -i ".$cuffdir."transcripts.gtf  -t ".$tmap." -c ".$cuffdir.$cuffcompare.".combined.gtf  -to ".$tmap_out."  -go $combined_out";

print "Selecting Cufflinks predictions to trim.\n";
system( $prepcmd );

##Run Trim_Cufflinks_GTF.pl 
my $trimcmd = "Trim_Cufflinks_GTF.pl -t ". $tmap_out." -g ".$combined_out." -c ".$coverage." -f ".$cutoff;

print"Trimming Cufflinks annotations. Water level is set at $cutoff of expected per-base coverage for gene.\n";
system($trimcmd);

##Run Get-TSS.pl
my $tsscmd = "Get-TSS.pl -t ".$tmap_out." -u ".$upstream." -d ".$downstream." -o ".$predictions_out;

print "Getting transcription start sites. Reporting $upstream bases upstream and $downstream bases downstream.\n";
system( $tsscmd );

##Run GetDistances.pl
my $distcmd = "GetDistances.pl -a ".$anno." -p ".$predictions_out." -o ".$distances_out;
print "Calculating distances between predicted and nearest annotated TSS.\n";

system( $distcmd );

print "Producing Distance and coverage cutoff plots";
system ( "TSSCutoffs.R" );

exit;
