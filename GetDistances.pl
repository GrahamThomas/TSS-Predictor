#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case no_auto_abbrev);

my ($annos,$preds,$outfile) = ("","tss.prediction.bed","distances.txt");
GetOptions (
  "a|anno:s" => \$annos,
  "p|predictions:s" => \$preds,
  "o|output:s" => \$outfile
);

##Read in the predictions, store TSS and annotation. Remember the order of the input TSS predictions (output file in the same order...)
my %prediction;
my @order;

open( PRED, $preds );

while(<PRED>){
	chomp;
	my @line = split( /\t/ );
	push( @order, $line[3] );
	if ($line[5] eq '+'){
		$prediction{$line[3]}{'PredStart'} = $line[1];
		$prediction{$line[3]}{'Chr'} = $line[0];
		$prediction{$line[3]}{'Gid'} = $line[3];
	} else {
		$prediction{$line[3]}{'PredStart'} = $line[2];
		$prediction{$line[3]}{'Chr'} = $line[0];
		$prediction{$line[3]}{'Gid'} = $line[3];
	}
}
close PRED;


##Go through the annotations and find nearest TSS to each gene.

open( ANNO, $annos );

my %what;

my $start;
while(<ANNO>){
	$_ =~ m/^.*?\t.*?\t(.*?)\t/;
	next unless $1 eq 'exon';
	$_ =~ m/gene_id \"(.*?)\".*?transcript_id \"(.*?)\".*?exon_number \"(\d+?)\"/;
	next unless $3 == 1; ##Only consider exon if is first exon in transcript.
	my $gid = $1;
	my $tid = $2;
	if ( defined $prediction{$gid} ){
		chomp;
		my @line = split(/\t/);
		if ( $line[6] eq '+'){ ##Get TSS from annotation
			$start = $line[3];
		} elsif ($line[6] eq '-') {
			$start = $line[4];
		}
		my $dist = $prediction{$gid}{'PredStart'} - $start;
		if ( not defined $prediction{$gid}{'Distance'} ){ ##Add nearest transcript, distance and TSS of nearest transcript information.
			$prediction{$gid}{'AnnoStart'} = $start;
			$prediction{$gid}{'Distance'} = $dist;
			$prediction{$gid}{'NearestTranscript'} = $tid;
		} elsif ( abs($prediction{$gid}{'Distance'}) > abs($dist) ) {
			$prediction{$gid}{'AnnoStart'} = $start;
			$prediction{$gid}{'Distance'} = $dist;
			$prediction{$gid}{'NearestTranscript'} = $tid;
		} else {
			next;
		}
	}
}
close ANNO;

##Print the output.
open( OUT, ">$outfile");

print OUT "Gene\tChr\tPrediction_start\tNearest_transcript\tNearest_transcript_start\tDistance\n";
foreach my $gene( @order ){
	my $chr = $prediction{$gene}{'Chr'};
	my $predst = $prediction{$gene}{'PredStart'};
	my $tran = $prediction{$gene}{'NearestTranscript'};
	my $transt = $prediction{$gene}{'AnnoStart'};
	my $dist = $prediction{$gene}{'Distance'};
	print OUT "$gene\t$chr\t$predst\t$tran\t$transt\t$dist\n";
}
close OUT;

exit;
