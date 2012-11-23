#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case no_auto_abbrev);

my ($anno,$isotigs,$tmap,$combinedgtf,$tmapout,$gtfout,$chr) = ("","transcripts.gtf","","","trim.gtf.tmap","trim.combined.gtf","");#These are the primary (ensembl etc) annotations, cufflinks' isotigs (transcripts.gtf), cufflinks' tmap file, cufflink's combined.gtf file, selected contigs from tmap output, trimmed combined.gtf output  
GetOptions (
  "a|anno:s" => \$anno,
  "i|isotigs:s" => \$isotigs,
  "t|tmap:s" => \$tmap,
  "c|combinedgtf:s" => \$combinedgtf,
  "to|tmapout:s" => \$tmapout,
  "go|gtfout:s" => \$gtfout,
  "ch|chromosome:s" =>\$chr,
);

###Parse the initial GTF annotation file (ensembl, etc) and links transcript ids to gene ids and gets the strand information.
my %tran2gene;
my %annostrand;

open( ANNOS, $anno );
while( <ANNOS> ){
	my @gtfline = split(/\t/, $_);
	if( $gtfline[8] =~ m/gene_id \"(.*?)\".*?transcript_id \"(.*?)\"/ ){
		$tran2gene{$2} = $1;
		$annostrand{$1} = $gtfline[6];
	} else {
		print "Could not successfully parse $anno. gene_id and transcript_id not in field 9. Are you sure this is a GTF file?\n";
		exit;
	}
}
close ANNOS;


###Parse the transcripts.gtf file, assign isotigs to strand Ignore anything not on $chr, if defined.
my %isotigstrand;

open( ISO, $isotigs );
NEXTISO: while( <ISO> ){
	my @gtfline = split(/\t/, $_);
#	if (defined $chr){
		next NEXTISO unless ( $chr eq $gtfline[0] || $chr eq "");
#	}
	$gtfline[8] =~ m/transcript_id \"(.*?)\"/;
	$isotigstrand{$1} = $gtfline[6]; 
}
close ISO;


###Read the transcripts.gtf.tmap file. Identify the longest cufflink's isotig for each gene. Get total gene length and coverage * length for each gene. Consider only transcripts ids present in %isotigstrand, this is $chr filter for tmap file.
my %longtransfrag;
my %tmap;
my %refgene;
my %cuffisodata;
my $refgene;

open( TMAP, $tmap );
my $tmap_head = <TMAP>; #Skip the header line
NEXTTMAP: while( <TMAP> ){
	my @line = split( /\t/, $_);
	next if( $line[2] eq 'u' || $line[2] eq 'r' ); ## Skip isotigs not assigned to genes
	next NEXTTMAP unless defined $isotigstrand{ $line[4] };
	$tmap{$line[4]} = $_; 
	$refgene = $tran2gene{ $line[1] };
	if( $annostrand{ $refgene } eq $isotigstrand{ $line[4] } ){ ###Strand filter.
			$cuffisodata{ $line[4] }{'refgene'} = $refgene; ###Get isotig data for all isotigs annotated to same strand as parent reference gene
			$cuffisodata{ $line[4] }{'length'} += $line[10];
			$cuffisodata{ $line[4] }{'coverage'} += $line[9];
		if( not defined $longtransfrag{ $refgene } ){ ###Identify the longest transfrag relating to each reference gene. (tmap file deals with transcripts, not genes!)
			$longtransfrag{$refgene}{'isotig'} = $line[4];
			$longtransfrag{$refgene}{'length'} = $line[10];
		} elsif ( $longtransfrag{ $refgene }{'length'} < $line[10] ) {
			$longtransfrag{$refgene}{'isotig'} = $line[4];
			$longtransfrag{$refgene}{'length'} = $line[10];
		}
	}
}
close( TMAP );

###Calculate expected coverage for entire gene as follows: (Iso1 coverage / longest_iso_length) * Iso1 length + (Iso2 coverage / longest_iso_length) * Iso2 length ... + (IsoN coverage / IsoN length) * IsoN length. Note - this is done at the level of the reference gene NOT the cufflinks gene.
##Also make a lookup table for the longest isotig for each refgene. 
my %normexp;
my %lookup;

foreach my $isotig( sort keys %cuffisodata ){
	$refgene = $cuffisodata{$isotig}{'refgene'};
	my $length = $cuffisodata{$isotig}{'length'};
	my $coverage = $cuffisodata{$isotig}{'coverage'};
	my $longestlength = $longtransfrag{$refgene}{'length'};
	$normexp{ $refgene } += ($coverage /$longestlength) * $length;
	$lookup{ $longtransfrag{$refgene}{'isotig'} } = '';
}


###Print a new tmap file with only longest isotigs for each reference gene and amended coverage estimates

open(OUT1, ">$tmapout");
my @tmap_out = split( /\t/, $tmap_head);
$tmap_out[9] = 'gene.cov';
print OUT1 join ( "\t", @tmap_out );

foreach my $outisotig( sort keys %longtransfrag ){
	my $isotig = $longtransfrag{ $outisotig }{ 'isotig' };
	my @out = split(/\t/, $tmap{ $isotig });
	$refgene = $tran2gene{ $out[1] };
	$out[1] = $refgene;
	my $normexp = $normexp{ $refgene };
	$out[9] = $normexp;
	print OUT1 join("\t", @out);
}
close OUT1;

open( ISO, "<$combinedgtf" );
open( OUT2, ">$gtfout");

while( <ISO> ){
	my $line = $_;
	$_ =~ m/oId \"(.*?)\"/;
	if (exists $lookup{ $1 } ){
		print OUT2 $line;
	}
}
close OUT2;


exit;
