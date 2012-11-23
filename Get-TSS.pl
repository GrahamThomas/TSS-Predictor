#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case no_auto_abbrev);

my ($gtf,$trimtmap,$upstream,$downstream, $outfile) = ("isotigs.trimmed.gtf","trim.gtf.tmap",0,1,"TSS.predictions.bed");
GetOptions (
  "g|gtf:s" => \$gtf,
  "t|trim-gtf-tmap:s" => \$trimtmap,
  "u|upstream:s" => \$upstream,
  "d|upstream:s" => \$downstream,
  "o|output:s" => \$outfile,
);

open( GTF, $gtf );

my %fgene;
my %rgene;
my %chr;
my ($gid, $tss, $chr, $strand, $start, $end);

my $count = 0;
while(<GTF>){
	my @line = split (/\t/, $_);
	$_ =~ m/^chr(.*?)\t.*?\t.*?\t(\d+?)\t(\d+?)\t\.\t(.*?)\t\.\t.*?transcript_id\s\"(.*?)\".*?contain_tss\s\"(.*?)\"/;
	$gid = $5;
	$tss = $6;
	$chr = $1;
	$strand = $4;
	$start = $2;
	$end = $3;	
	if ($strand eq "+"){
		if ($tss eq "yes"){
			$fgene{$chr}{$gid} = $start; 
		}
	} elsif ( $strand eq "-" ){
		if ($tss eq "yes"){
			$rgene{$chr}{$gid} = $end;
		}
	}
}


open( TMAP, $trimtmap );
my %cuff2ref;
while(<TMAP>){
	my @line = split( /\t/, $_ );
	$cuff2ref{ $line[4] } = $line[1];
}
close TMAP;


open (OUT,">$outfile" );
#Print out locations of the promoter sequences for genes on the + strand
foreach $chr( sort { $fgene{$a} cmp $fgene{$b} } keys ( %fgene )  ){
	foreach my $pos( sort{ $fgene{$chr}{$a} <=> $fgene{$chr}{$b} } keys ( %{$fgene{$chr}} ) ){
		my $us = $fgene{$chr}{$pos} - $upstream;
		my $ds = $fgene{$chr}{$pos} + $downstream;
 		print OUT "chr".$chr."\t".$us."\t".$ds."\t".$cuff2ref{ $pos }."\t\.\t+\n"; ###Print as BED
	}
}

#Print out locations of the promoter sequences for genes on the - strand
foreach $chr( sort { $rgene{$a} cmp $rgene{$b} } keys ( %rgene )  ){
	foreach my $pos( sort{ $rgene{$chr}{$a} <=> $rgene{$chr}{$b} } keys ( %{$rgene{$chr}} ) ){
		my $us = $rgene{$chr}{$pos} + $upstream;
		my $ds = $rgene{$chr}{$pos} - $downstream;
 		print OUT "chr".$chr."\t".$ds."\t".$us."\t".$cuff2ref{ $pos }."\t\.\t-\n"; ###Print as BED, downstream and upstream the opposite way around due to opposite strand!
	}
}

close OUT;
exit;
