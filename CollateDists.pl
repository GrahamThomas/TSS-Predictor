#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my %dists;
my $outdir = shift(@ARGV);
my @order;

###Iterate through the list of distances.txt files and parse the cutoff used as well as the distance or each gene.
foreach my $file(@ARGV){
	$file =~ m/^$outdir\/(.*?).distances.txt/;##Get the cutoff used
	my $cutoff = $1;
	push ( @order, $cutoff );
	open( FILE, "<$file");
	<FILE>;#Skip the header line
	while(<FILE>){
		my @line = split(/\t/);
		my $gene = $line[0];#Get the reference gene ID
		my $dist = $line[5];#Get the distance from prediction to nearest reference gene TSS
		chomp $dist;
		$dists{ $gene }{ $cutoff } = $dist;
	}
	close FILE;
}

open( OUT, ">$outdir/CutoffDistances.txt" );
my $header = "Gene_id\t";
foreach my $cutoff( @order ){
	$header .= "$cutoff\t";	
} 

chop $header;
print OUT "$header\n";

foreach my $gene( sort { $a cmp $b } keys %dists ){
	print OUT "$gene";
	foreach my $cutoff( @order ){
		if( defined $dists{$gene}{$cutoff} ){		
			print OUT "\t$dists{$gene}{$cutoff}";
		} else {
			print OUT "\tNA";
		}
	}
	print OUT "\n";
} 

exit;
