#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case no_auto_abbrev);
use Data::Dumper;

my ($distance, $expression, $outfile, $refout ) = (undef,undef,"tss.selected","0");
GetOptions (
  "d|distance:s" => \$distance,
  "e|expr:s" => \$expression,
  "o|output:s" => \$outfile,
  "r|refout" => \$refout,###Make this a binary tag, rather than string.
);

my $tmap = "trim.gtf.tmap";
my $predictions = "tss.prediction.bed";
my $distfile = "distances.txt";

my $usage = "SelectTSS.pl [options -d -e -r] -o <output prefix>
Options:-
Short	Long	Default					Description
-e	expr		[Optional]			Per-base coverage (expression) cutoff.
-d	distance	[Optional]			Distance between predicted and nearest annotated TSS to report
-r	refout		[Optional]			Report nearest reference TSS rather than the predicted TSS
-o	output		[tss.selected]			Output file prefix
";

die "$usage\n" unless $predictions && $distfile && $tmap && ( defined($distance)||defined($expression)||($refout==1));

if( defined $expression ){
	die "$usage\n\nExpression cutoff '-e $expression' not numeric.\n" unless $expression =~ m/^[0-9|\.]+$/;
}

if( defined $distance ){
	die "$usage\n\nDistance cutoff '-d $distance' not numeric.\n" unless $distance =~ m/^[0-9|\.]+$/;
}

##Create a blacklist. The all genes that fail any filter are added.
my %blacklist;
my %distance;
my @line;

###Parse accepted genes based on expression
if( defined $expression ){
	die "Error! Could not open $tmap file.\n" unless open(TMAP, $tmap);
	<TMAP>;
	while(<TMAP>){
		@line = split("\t");
		if( $line[9] < $expression ){
			$blacklist{$line[1]} = 0;
		}
	}
}
close TMAP;


###If either distance or -r flag provided then parse the distances.txt file.
my $genedist;

if( (defined $distance) || ($refout == 1) ){
	die "Error! Could not open $distfile file.\n" unless open( DISTS, $distfile );
	<DISTS>;
	while(<DISTS>){
		@line = split("\t");
		if( defined $distance ){
			$line[5] =~ m/(\d+?)\n/;
			if( $1 > $distance ){
				$blacklist{$line[0]}= 0;
			}
		}
		if( $refout == 1 ){
			$distance{$line[0]} = $line[5];
		} 
	}
}
close DISTS;

##Parse through distances file again, print only lines for genes not on blacklist
open( DISTOUT, ">$outfile\.distances.txt" );
open( DISTS, $distfile );
while(<DISTS>){
	@line = split("\t");
	if( defined $blacklist{$line[0]} ){
		next;
	}	else {
		print DISTOUT;
	}
}
close DISTS;
close DISTOUT;

##Now do the same for the TSS predictions, amend the TSS position by the offset in %distance if necessary.
open( TSSOUT, ">$outfile\.bed" );
die "Error! Could not open $predictions file.\n" unless open( TSS, "<$predictions" );
LINE:while(<TSS>){
	@line = split("\t");
	if( defined $blacklist{$line[3]} ){
		next LINE;
	} else {
		if( $refout == 1 ){###Amend TSS coordinates to reference if aplicable.
			$line[1] = $line[1] - $distance{$line[3]};
			$line[2] = $line[2] - $distance{$line[3]};
			}
	my $outline = join( "\t", @line );
	print TSSOUT "$outline";
	}
}

exit;
