#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

use Getopt::Long qw(:config pass_through no_ignore_case no_auto_abbrev);

my ($ref,$ensembl,$out) = ("","","");
GetOptions (
  "r|ref:s" => \$ref,
  "e|ensembl:s" => \$ensembl,
  "o|out:s" => \$out,
);

my $usage = "Ensembl2UCSC.pl -r <reference> -e <Ensembl GTF> -o <Output file>
Options:-
Short	Long			Description
-r	reference		Multi-FASTA file containing the reference genome
-e	ensembl			Ensembl GTF file
-o	output			Output file	
";

die "$usage\n" unless $ref && $ensembl && $out;

my @reference;

open(REF, $ref);
while(<REF>){
	next unless $_ =~ m/^>(chr.*?)\n/;
	push (@reference, $1);
}


my @annotation;
my $line;

open( ENS, $ensembl );
while(<ENS>){
	$line = $_;
	$line =~ s/^MT/M/;
	$line =~ s/(.*?)\s/chr$1\t/;
	push(@annotation, $line);
}

#Iterate successively through the output order file and the annotation hash and print the line if the first column matches the output order. This is inefficient as it requires going throung the annotation hash multiple times however is a lazy way of ensuring that annotations are kept in their chromosomal order.

my $chr;

open( OUT, ">$out" );
foreach my $outorder( @reference ){
	foreach my $annoline( @annotation ){
		$annoline =~ m/(.*?)\t/;
		$chr = $1;
		if ( $outorder eq $chr){
			print OUT $annoline;
		}
	}
}
close OUT;

exit;

