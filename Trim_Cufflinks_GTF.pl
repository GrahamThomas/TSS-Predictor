#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case no_auto_abbrev);

my ($interest,$combined,$coverage,$cutoff,$outfile,$chrom) = ("trim.gtf.tmap","trim.combined.gtf","","0.1","isotigs.trimmed.gtf",""); 
GetOptions (
  "t|tmap:s" => \$interest,
  "g|combined-gtf:s" => \$combined,
  "c|coverage:s" => \$coverage,
  "f|cutoff_fraction:s" => \$cutoff,
  "o|output:s" => \$outfile,
  "ch|chromosome" => \$chrom,
);

my %gene_info;
my %exon_info;
my %coverage;
my %is_interest;

##Read in the isoforms of interest and the expected coverage

open( INTEREST, $interest );
<INTEREST>; ###Skip the header line
while(<INTEREST>){
	my @line = split(/\t/, $_);
	$is_interest{$line[4]} = $line[9]; ###Get the transcript id of interest (key) and expected coverage (value)
}
close INTEREST;

##Read in the gene information of the genes of interest from the combined.gtf file

my ($chr, $start, $end, $strand, $gene_id, $transcript_id, $exon_number, $gene_name, $oid, $nearest, $class_code, $tss_id); 
my $exon_counter = 0;
my %exon_lookup; ##hash table to contain a link between the unique exon id and the transcript it belongs to

open(GTF, $combined);
while(<GTF>){
	my @line = split(/\t/, $_);
	$line[8] =~ m/gene_id\s\"(.*?)\";\stranscript_id\s\"(.*?)\";\sexon_number\s\"(.*?)\";\sgene_name\s\"(.*?)\";\soId\s\"(.*?)\";\snearest_ref\s\"(.*?)\";\sclass_code\s\"(.*?)\";\stss_id\s\"(.*?)\".*?\n/;
	$chr=$line[0];
	$start = $line[3];
	$end = $line[4];
	$strand = $line[6];
	$gene_id = $1;
	$transcript_id = $2;
	$exon_number = $3;
	$gene_name = $4;
	$oid = $5;
	$nearest = $6;
	$class_code = $7;
	$tss_id = $8;
	if (defined $is_interest{ $oid } ){
		$gene_info{$chr}{$oid}{'strand'} = $strand;
		$gene_info{$chr}{$oid}{'gene_id'} = $gene_id;
		$gene_info{$chr}{$oid}{'gene_name'} = $gene_name;
		$gene_info{$chr}{$oid}{'strand'} = $strand;
		$gene_info{$chr}{$oid}{'oId'} = $oid;
		$gene_info{$chr}{$oid}{'transcript_id'} = $transcript_id;
		$gene_info{$chr}{$oid}{'nearest_ref'} = $nearest;
		$gene_info{$chr}{$oid}{'class_code'} = $class_code;
		$gene_info{$chr}{$oid}{'tss_id'} = $tss_id;
		$gene_info{$chr}{$oid}{'coverage'} = $is_interest{ $oid };
		
		my $exon_id = "$exon_counter";
		$exon_info{$chr}{$exon_id}{'exon_id'} = $exon_counter;
		$exon_info{$chr}{$exon_id}{'oId'} = $oid;
		$exon_info{$chr}{$exon_id}{'exon_number'} = $exon_number;
		$exon_info{$chr}{$exon_id}{'start'} = $start;
		$exon_info{$chr}{$exon_id}{'end'} = $end;
		$exon_info{$chr}{$exon_id}{'coverage'} = $is_interest{ $oid };
		$exon_info{$chr}{$exon_id}{'exon_code_f'} = 2; ##Add this variable. Later amend this to state whether whole exon is kept (2 = default), trim exon on f strand = 1, disregard exon entirely = 0.
		$exon_info{$chr}{$exon_id}{'exon_code_r'} = 2; ##Add this variable. Later amend this to state whether whole exon is kept (2 = default), trim exon on r strand = 1, disregard exon entirely = 0.		
		$exon_info{$chr}{$exon_id}{'contain_tss'} = "no"; ##By default all exons have tss = no. When those containing tss's are identified later this is amended to yes.
		$exon_counter++;		
		
		$exon_lookup{ $chr }{ $exon_id } = $oid;
		}
	
}
close GTF;

###Read in the coverage information

open( COV, $coverage );
while(<COV>){
	my @line = split (/\t/, $_);
	chomp $line[2];
	$coverage{ $line[0] }{ $line[1] } = $line[2]; ###Coverge{ Chr }{ Pos } = Coverage
}
close COV;

##Iterate over the genes (indicies of @exon_start) and ask whether coverage at each base in an exon is $cutoff proportion of the expected coverage. 
##Once this condition is met add the variable $gene_info{$chr}{$transcript_id}{'gene_bound_f'} containing the start pos of the gene, skip all exons belonging to genes once flag is assigned.

my ( $exon_start, $exon_end, $tid, $exp_cov, $loc, $exon_id );
my $lastloc =0;

foreach my $chr( sort { $exon_lookup{$a} cmp $exon_lookup{$b} } keys ( %exon_lookup ) ){
	EXON: 
	foreach my $exon_sort( sort { $a <=> $b } keys ( %{$exon_lookup{$chr}} ) ){
		$exon_start = $exon_info{$chr}{$exon_sort}{'start'};
		$exon_end = $exon_info{$chr}{$exon_sort}{'end'};
		$tid = $exon_info{$chr}{$exon_sort}{'oId'};
		$exp_cov = $exon_info{$chr}{$exon_sort}{'coverage'};
		my $istrim =0;
		if ( defined $gene_info{$chr}{$tid}{'gene_bound_f'} ){ ##Skip all exons for gene once end position has been amended
			next EXON;
		}
		LOC:
		for $loc( reverse( $exon_start .. $exon_end ) ){
			if (not defined $coverage{$chr}{$loc}){$coverage{$chr}{$loc} = 0;} ##pileup output doesn't report bases with 0 coverage. add these values here.
			if ($coverage{$chr}{$loc} > ($exp_cov * $cutoff) ){
				$istrim = 1;
				$lastloc = $loc;
				if( $exon_info{$chr}{$exon_sort}{'start'} == $loc ){
					$exon_info{$chr}{$exon_sort}{'exon_code_f'} = 1;
					$gene_info{$chr}{$tid}{'gene_bound_f'} = $loc;
					next EXON;
				}
				next LOC;
			} 
			elsif ($coverage{$chr}{$loc} < ($exp_cov * $cutoff) ){
				if ( $istrim == 0 ){
					for my $loc2( reverse( $exon_start .. $exon_end ) ){ ###If the expression at the end of a first exon is < the cutoff then go 20% into the exon. If at no point in this window does the coverage pop above the cutoff then disregard the exon, otherwise include it in the analysis.
						if ( $coverage{$chr}{$loc2} > ($exp_cov * $cutoff) ){
							next LOC;													
						} else {
							$exon_info{$chr}{$exon_sort}{'exon_code_f'} = 0;
							next LOC;						
						}
					}
				}
				if ($lastloc > ( $exp_cov * $cutoff ) ){
					$exon_info{$chr}{$exon_sort}{'exon_code_f'} = 1;
					$gene_info{$chr}{$tid}{'gene_bound_f'} = $lastloc;
					next EXON;
				}
			}
#			if ( ($istrim == 1) && (not defined $exon_info{$chr}{$exon_sort}{'exon_code_f'} ) ){###This is supposed to ensure that I don't walk off the end of a gene without reporting it. It doesn't work as exon_code_f is always defined. Fix is inserted into loop if ($coverage{$chr}{$loc} > ($exp_cov * $cutoff) ){... Assign gene_bound_r if we walk into the last base of the cufflinks annotation
#				$exon_info{$chr}{$exon_sort}{'exon_code_f'} = 1;
#				$gene_info{$chr}{$tid}{'gene_bound_f'} = $lastloc;
#				next EXON;
#			}
		}
	}
}


###Now, do the same in the reverse order. This is all working well!

foreach my $chr( sort { $exon_lookup{$a} cmp $exon_lookup{$b} } keys ( %exon_lookup ) ){
	EXON: 
	foreach my $exon_sort( sort { $b <=> $a } keys ( %{$exon_lookup{$chr}} ) ){
		$exon_start = $exon_info{$chr}{$exon_sort}{'start'};
		$exon_end = $exon_info{$chr}{$exon_sort}{'end'};
		$tid = $exon_info{$chr}{$exon_sort}{'oId'};
		$exp_cov = $exon_info{$chr}{$exon_sort}{'coverage'};
		my $istrim =0;
		if ( defined $gene_info{$chr}{$tid}{'gene_bound_r'} ){ ##Skip all exons for gene once end position has been amended
			next EXON;
		}
		LOC:
		for $loc( $exon_start .. $exon_end ){
			if (not defined $coverage{$chr}{$loc}){$coverage{$chr}{$loc} = 0;} ##pileup output doesn't report bases with 0 coverage. add these values here.
			if ($coverage{$chr}{$loc} > ($exp_cov * $cutoff) ){
				$istrim = 1;
				$lastloc = $loc;
				if( $exon_info{$chr}{$exon_sort}{'end'} == $loc ){
					$exon_info{$chr}{$exon_sort}{'exon_code_r'} = 1;
					$gene_info{$chr}{$tid}{'gene_bound_r'} = $loc;
					next EXON;
				}
				next LOC;
			} 
			elsif ($coverage{$chr}{$loc} < ($exp_cov * $cutoff) ){
				if ( $istrim == 0 ){
					for my $loc2( $exon_start .. $exon_end ){ ###If the expression at the end of a first exon is < the cutoff then go 20% into the exon. If at no point in this window does the coverage pop above the cutoff then disregard the exon, otherwise include it in the analysis.
						if ( $coverage{$chr}{$loc2} > ($exp_cov * $cutoff) ){
							next LOC;													
						} else {
							$exon_info{$chr}{$exon_sort}{'exon_code_r'} = 0;
							next LOC;						
						}
					}
				}
				if ($lastloc > ( $exp_cov * $cutoff ) ){
					$exon_info{$chr}{$exon_sort}{'exon_code_r'} = 1;
					$gene_info{$chr}{$tid}{'gene_bound_r'} = $lastloc;
					next EXON;
				}
			}
#			if ( ($istrim == 1) && (not defined $exon_info{$chr}{$exon_sort}{'exon_code_r'} ) ){###This is supposed to ensure that I don't walk off the end of a gene without reporting it. It doesn't work as exon_code_f is always defined. Fix is inserted into loop if ($coverage{$chr}{$loc} > ($exp_cov * $cutoff) ){... Assign gene_bound_r if we walk into the last base of the cufflinks annotation
#				$exon_info{$chr}{$exon_sort}{'exon_code_r'} = 1;
#				$gene_info{$chr}{$tid}{'gene_bound_r'} = $lastloc;
#				next EXON;
#			}
		}
	}
}

####Amend the start and end positions of the exons where appropriate. Print exons with 'exon_code's as follows: 0-do not print, 1-amend position (gene_bound_f for exon_code_f =1, gene_bound_r for exon_code_r =1), 3-print as is. 
#Also add the flag "contains_tss" = yes to exons containing exon_code 1/2 and "contain_tss" = no when exon_code = 3. 

foreach my $chr( sort { $exon_lookup{$a} cmp $exon_lookup{$b} } keys ( %exon_lookup ) ){
	foreach my $exon_sort( sort { $a <=> $b } keys ( %{$exon_lookup{$chr}} ) ){
		###Amend the start and end positions of the exon (if necessary)
		my $tid = $exon_info{$chr}{$exon_sort}{'oId'};
		if ( $exon_info{$chr}{$exon_sort}{'exon_code_f'} == 0 || $exon_info{$chr}{$exon_sort}{'exon_code_r'} == 0 ){
			delete ($exon_info{$chr}{$exon_sort}); ###Delete this exon record entirely 
		} elsif ( $exon_info{$chr}{$exon_sort}{'exon_code_f'} == 1 ){
			$exon_info{$chr}{$exon_sort}{'start'} = $gene_info{$chr}{$tid}{'gene_bound_f'};
			if ( $gene_info{$chr}{ $exon_info{$chr}{$exon_sort}{'oId'} }{'strand'} eq "+" ){
				$exon_info{$chr}{$exon_sort}{'contain_tss'} = "yes";			
			}
		} elsif ( $exon_info{$chr}{$exon_sort}{'exon_code_r'} == 1 ){
			$exon_info{$chr}{$exon_sort}{'end'} = $gene_info{$chr}{$tid}{'gene_bound_r'};
			if ( $gene_info{$chr}{ $exon_info{$chr}{$exon_sort}{'oId'} }{'strand'} eq "-" ){
					$exon_info{$chr}{$exon_sort}{'contain_tss'} = "yes";
			}		
		} 
	}
}


###Finally reconstitute the GTF file.
open( OUT, ">$outfile");
foreach my $chr( sort { $exon_lookup{$a} cmp $exon_lookup{$b} } keys ( %exon_lookup ) ){
	foreach my $exon_sort( sort { $a <=> $b } keys ( %{$exon_lookup{$chr}} ) ){
		if (not defined $exon_info{$chr}{$exon_sort}{'oId'}){
			next;
		}
		my $oid = $exon_info{$chr}{$exon_sort}{'oId'};
		$start = $exon_info{$chr}{$exon_sort}{'start'};
		$end = $exon_info{$chr}{$exon_sort}{'end'};
		$exon_number = $exon_info{$chr}{$exon_sort}{'exon_number'};
		my $tss = $exon_info{$chr}{$exon_sort}{'contain_tss'};
		$strand = $gene_info{$chr}{$oid}{'strand'};
		$gene_name = $gene_info{$chr}{$oid}{'gene_name'};
		$nearest = $gene_info{$chr}{$oid}{'nearest_ref'};
		$coverage = $gene_info{$chr}{$oid}{'coverage'};
		print OUT "$chr\tCufflinks_trimmed\texon\t$start\t$end\t.\t$strand\t.\tgene_id \"$oid\"; transcript_id \"$oid\"; exon_number \"$exon_number\"; gene_name \"$gene_name\"; nearest_ref \"$nearest\"; coverage \"$coverage\"; contain_tss \"$tss\";\n";
	}
}

close OUT;
exit;
