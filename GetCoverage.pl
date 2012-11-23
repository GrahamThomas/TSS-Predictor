#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case no_auto_abbrev);

my ($combined, $bam, $genome, $coverage ) = ("","","","");
GetOptions (
  "c|combined:s" => \$combined,
  "b|bam:s" => \$bam,
  "g|genome:s" => \$genome,
  "o|output:s" => \$coverage,
);

my $usage = "GetCoverage.pl -c <cufflinks.combined.gtf> -b <mapped_reads.bam> -g <.genome file> -o <output>
Options:-
Short	Long		Default			Description
-c	combined	[]			This is the Cuffcompare combined.gtf file
-b	bam		[]			A sorted BAM file
-g	genome		[]			A .genome file as provided by BEDTools
-o	output		[]			output .coverage file
";

die "$usage\n" unless $combined && $bam && $genome && $coverage;


##Put some dies in if the file handles aren't there.


my @args = ( "bedtools --version" );
system( @args ) == 0 or die "could not find BEDTools. Are you sure it is properly installed?";
print "BEDTools found\n";

##Create a merged GTF file
my $merged = $combined;
$merged =~ s/^(.*?combined)(\.gtf)$/$1\.merged.bed/;

my $mergecmd = "sortBed -i $combined | mergeBed -i stdin > $merged";
system ( $mergecmd );
if ( $? == -1 ){
	print "could not merge the combined.gtf file. Perhaps BedTools isn't properly installed?";
}

print "Calculating per-base coverage, this may take some time.\n";

my $coveragecmd = "intersectBed -wa -abam $bam -b $merged | genomeCoverageBed -d -split -ibam stdin -g $genome |". ' perl -ne \'m/^(\w+?)\t(\d+?)\t(\d+?)\n/; $end=$2+1; print"$1\t$2\t$end\t$3\n"\' - |' ."intersectBed -wa -a stdin -b $merged | ".'perl -ne \'m/^(\w+?)\t(\d+?)\t\d+?\t(\d+?)\n/; print"$1\t$2\t$3\n"\''." > $coverage";
system ($coveragecmd);

exit;

#GetCoverage.pl -c chr1.multi.ensGene58.cuffcompare.combined.gtf -b ~/Desktop/MacrophageTSS/TopHat/Macrophage.chr1.bam -g ~/Programs/bedtools/BEDTools-Version-2.16.2/genomes/mouse.mm9.genome -o test.coverage