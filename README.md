TSS-Predictor
=============

Version 1.0 
by Graham Thomas


-----------------


Contents
=======

###1.0. Introduction###

###2.0. Getting Started###
	2.1. Dependencies
	2.2. Installation
	
###3.0. Using TSS-Predictor###
	3.1. Preparing your data
	3.2. Mapping reads to the reference genome
	3.3. Running Cufflinks and Cuffcompare
	3.4. Generating a .coverage file with GetCoverage.pl
	3.5. Using TSS-Optimizer to identify optimal coverage cutoffs
	3.6. Running TSS-Predictor.pl
	3.7. Post-filtering for sensible results
	3.8. Seecting a high confidence TSS set using SelectTSS.pl

###4.0. Output files and formats###
	4.1. Coverage file
	4.2. trim.combined.gtf
	4.3. trim.gtf.tmap
	4.4. isotigs.trimmed.gtf
	4.5. tss.prediction.bed
	4.6. distances.txt
	4.7. tss.selected.bed

###5.0. Overview of the algorithm###
	5.1. Prep4Trim.pl
	5.2. Trim_Cufflinks_GTF.pl
	5.3. Get-TSS.pl
	5.4. GetDistances.pl
	5.5. SelectTSS.pl
	5.6. CollateDistances.pl


-----------------


1.0 Introduction
================

TSS-Predictor has been developed to identify promoter usage for genes in RNA-Seq experiments. It is designed to be used when both a high quality reference genome and annotation set are available, for example human or mouse. The output of TSS-Predictor is a BED file containing predicted transcription start sites for each reference gene, these may then be used for downstream applications such as _cis_-regulatory analysis.

TSS-Predictor is written in perl and is available at GitHub <https://github.com/GrahamThomas/TSS-Predictor> and as an online tool at Galaxy (Not Yet!!!) and GeneProf (Not Yet!!!).


2.0. Getting Started
====================

2.1. Dependencies
-----------------

In order to use TSS-Predictor you must have working versions of the following in your $PATH:- 

1. A short read mapper compatible with Cufflinks. I use TopHat <http://tophat.cbcb.umd.edu/>.
2. Cufflinks <http://cufflinks.cbcb.umd.edu/>
3. BEDTools <http://code.google.com/p/bedtools/> and
4. R, with the libraries ggplot2 and reshape2 installed <http://cran.r-project.org/>

To follow this guide you will need to merge and sort some BAM files. For this you will also need SAMtools <http://samtools.sourceforge.net/>. 

2.2. Installation
-----------------

TSS-Predictor is implemented in Perl and R. To get up and running all you need to do is:-

1. Add TSS-Predictor/bin to your $PATH variable:-

	\$PATH=$PATH:/\<Path to installation directory\>/TSS-Predictor/bin

	

2. And give executable permission to All scripts in the directory.


	cd \<Path to installation directory\>/TSS-Predictor/bin

	chmod +x \*


-----------------

3.0. Using TSS-Predictor
========================

3.1. Preparing your data
------------------------

Below is a step-by-step guide on how to use TSS-Predictor beginning with raw illumina reads. In this example we use 51 base paired-end Illumina RNA-Seq data from a macrophage transcriptomics project (REF). Raw data and BAM files for this project are available from the SRA (ACCESSION... <http://www.ebi.ac.uk/ena/>). All of the intermediate files generated in the working example are available at our FTP site <link...>.


 	cd ~/TSS-Predictor-Example/RawData/
 	ls
 	4Ne_1_F.fastq  4Ne_3_R.fastq  4TG_3_F.fastq  BNe_2_R.fastq  BTG_2_F.fastq            
	4Ne_1_R.fastq  4TG_1_F.fastq  4TG_3_R.fastq  BNe_3_F.fastq  BTG_2_R.fastq  
	4Ne_2_F.fastq  4TG_1_R.fastq  BNe_1_F.fastq  BNe_3_R.fastq  BTG_3_F.fastq  
	4Ne_2_R.fastq  4TG_2_F.fastq  BNe_1_R.fastq  BTG_1_F.fastq  BTG_3_R.fastq
 	4Ne_3_F.fastq  4TG_2_R.fastq  BNe_2_F.fastq  BTG_1_R.fastq
 
 
A pre-built Bowtie index for mouse (mm9) was downloaded from <http://bowtie-bio.sourceforge.net/index.shtml>.

	cd ~/TSS-Predictor-Example/BowtieIndexes/
	wget ftp://ftp.cbcb.umd.edu/pub/data/bowtie_indexes/mm9.ebwt.zip
	unzip mm9.zip

Reconstitute a reference genome from the bowtie index.

	cd ~/TSS-Predictor-Example/BowtieIndexes/
	bowtie-inspect mm9 > mm9.fa


An Ensembl GTF file containing reference annotations was downloaded from the Ensembl FTP site <ftp://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/>. Ensembl annotations are preferable as, unlike UCSC annotations, a many-to-one relationship exists between gene 'gene\_id' and 'transcript\_id' fields of the GTF file. This relationship is required to combine per-transcript expression levels to per-gene expression levels.
	
	cd ~/TSS-Predictor-Example/Annotations/	
	wget ftp://ftp.ensembl.org/pub/release-67/gtf/\
	mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz 

Ensembl chromosomes need converting from NCBIM37 to mm9 style. This can be done with Ensembl2UCSC.pl. This script prints the annotations in the same order as the FASTA headers in the supplied genome and removes annotations belonging to the 'random' chromosomes. An alternative approach would be to map reads directly to the Ensembl genome build, thereby ommitting this step.

	Ensembl2UCSC.pl -r ../BowtieIndexes/mm9.fa -e Mus_musculus.NCBIM37.67.gtf \ 
	-o Mus_musculus.mm9.67.gtf


3.2. Mapping reads to the reference genome
------------------------------------------
 
Map reads to the reference genome using TopHat (version v2.0.0, Bowtie version 0.1.18.0 were used in this example). Multi-mapping reads are allowed by using the default parameters. Multi-mapped reads are important as they facilitate the assembly of transcripts containing repetitive sequence. In downstream analysis it may be wise to disregard transcripts containing a large proportion of multi-mapped reads as they may be unreliable. 

Transcription start site prediction is more effective with greater sequencing depth. My preferred strategy is to map reads for each lane individually, and then merge the alignments prior to transcript assembly with Cufflinks.

Map the reads with TopHat:-	

	cd ~/TSS-Predictor-Example/TopHat
	tophat -p 5 -o KNe_1_multi_out ../BowtieIndexes/mm9 \
	../RawData/4Ne_1_F.fastq ../RawData/4Ne_1_R.fastq
	...
	tophat -p 5 -o BTG_3_multi_out ../BowtieIndexes/mm9 \
	../RawData/BTG_3_F.fastq ../RawData/BTG_3_R.fastq

 Merge and sort the alignments. I use SAMtools <http://samtools.sourceforge.net/>:-
 
	cd ~/TSS-Predictor-Example/TopHat/

	###Reconstitute a SAM header from one of the TopHat alignments
	samtools view -H BNe_1_multi_out/accepted_hits.bam > mm9.header.sam


	samtools cat -h mm9.header.sam -o Macrophage.cat.bam \
	BNe_1_multi_out/accepted_hits.bam BNe_2_multi_out/accepted_hits.bam \
	BNe_3_multi_out/accepted_hits.bam KNe_1_multi_out/accepted_hits.bam \ 
	KNe_2_multi_out/accepted_hits.bam KNe_3_multi_out/accepted_hits.bam \
	BTG_1_multi_out/accepted_hits.bam BTG_2_multi_out/accepted_hits.bam \
	BTG_3_multi_out/accepted_hits.bam KTG_1_multi_out/accepted_hits.bam \
	KTG_2_multi_out/accepted_hits.bam KTG_3_multi_out/accepted_hits.bam

	samtools sort Macrophage.cat.bam Macrophage.sorted
	
3.3. Running Cufflinks and Cuffcompare
--------------------------------------

Assemble the mapped reads into transcribed fragments 'transfrags' using Cufflinks (here we have used version 2.0.2). In this example we have not provided reference annotations to guide transcript assembly. You may do so if you wish.

	cd ~/TSS-Predictor-Example
	cufflinks -p 8 -L Macrophage-Cufflinks -o Cufflinks TopHat/Macrophage.sorted.bam


Assign the transfrags to reference gene models using Cuffcompare

	cd ~/TSS-Predictor-Example/Cufflinks
	cuffcompare -r ../Annotations/Mus_musculus.mm9.67.sorted.gtf -M \
	-s ../mm9.fa -o mm9-EnsGene transcripts.gtf

 
3.4. Generating a .coverage file with GetCoverage.pl
----------------------------------------------------

GetCoverage.pl makes use of BEDTools to calculate the per-base read coverage at the genomic locations reported in the combined.gtf file produced by Cuffcompare. The script takes the following arguments:-

	GetCoverage.pl -c <cuffcompare.combined.gtf> -b <mapped_reads.bam> \
	-g <.genome file> -o <output>

+---------+--------------------------------------------+
| Options | Description                                |
+=========+============================================+
| -c      | Path to the cuffcompare combined.gtf file  |
+---------+--------------------------------------------+
| -b      | Path to the mapped reads BAM file          |
+---------+--------------------------------------------+
| -g      | Path to the .genome file                   |
+---------+--------------------------------------------+
| -o      | Name of the output .coverage file          |
+---------+--------------------------------------------+


The .genome file is a two column tab-delimited file stating each chromosome/scaffold name and it's size in base pairs. BEDTools is shipped with .genome files for mouse and human reference genomes. If you are working with a different species you will have to create one yourself, see the BEDTools manual (<http://code.google.com/p/bedtools/>) for more information. 

Making the .coverage file can take some time. We make the .coverage file for this example as follows:-

	cd ~/TSS-Predictor-Example/CoverageFiles
	cp <Path To BEDTools>/genomes/mouse.mm9.genome .
	
	GetCoverage.pl -c ../Cufflinks/mm9-EnsGene.combined.gtf \
	-b ../TopHat/Macrophage.sorted.bam -g mouse.mm9.genome \
	-o mm9-EnsGene.coverage


3.5. Using TSS-Optimizer.pl to identify optimal cutoffs
-------------------------------------------------------

TSS-Optimizer.pl is a wrapper script for the main TSS prediction algorithm. For a description of how the main algorithm works see section 5 or the supplementary information of (REF MY PAPER). TSS-Optimizer.pl takes one chromosome, defined by the user and predicts transcription start sites over a range of cutoff fractions. TSS-Optimizer.pl reports, for each gene at each cutoff fraction, the distance between the predicted TSS and nearest annotated TSS for the same gene. TSS-Optimizer also produces a graph showing the median absolute deviation (MAD) between the predicted and nearest annotated transcription start sites at each cutoff fraction. We choose the cutoff fraction that gives the lowest MAD for use with TSS-Predictor.

	TSS-Optimizer.pl -a <Annotations.gtf> -cd <Cufflinks_directory> \
	-cc <Cuffcompare_prefix> -t <tmap_file> -co <Coverage_file> -ch <chromosome>

+---------+---------------+--------------------------------------------------+
| Options | Default       | Description                                      |
+=========+===============+==================================================+
| -a      | []            | The main reference annotation file               |
|         |               | i.e. the one used with Cuffcompare               |
+---------+---------------+--------------------------------------------------+
| -cd     | []            | Path to the Cufflinks directory                  |
+---------+---------------+--------------------------------------------------+
| -cc     | []            | Prefix used when running Cuffcompare             |
+---------+---------------+--------------------------------------------------+
| -t      | []            | Path to the Cuffcompare transcripts.gtf.tmap file|
+---------+---------------+--------------------------------------------------+
| -co     | []            | Path to the .coverage file produced by           |
|         |               | GetCoverage.pl                                   |
+---------+---------------+--------------------------------------------------+
| -ch     | []            | Name of chromosome for optimization              |
+---------+---------------+--------------------------------------------------+
| -s      | [0.01]        | Coverage cutoff fraction to begin trimming       |
+---------+---------------+--------------------------------------------------+
| -e      | [0.20]        | Coverage cutoff fraction to end trimming         |
+---------+---------------+--------------------------------------------------+
| -i      | [0.01]        | Increment, or step size, between lower and upper |
|         |               | coverage cutoff fractions                        |
+---------+---------------+--------------------------------------------------+
| -o      | [./Optimizer] | Output directory (this will be suffixed with the |
|         |               | selected chromosome)                             |
+---------+---------------+--------------------------------------------------+

Running TSS-Optimizer for our data set is simple. We use the following command:-

	cd ~/TSS-Predictor-Example
	
	TSS-Optimizer.pl -a Annotations/Mus_musculus.mm9.67.gtf \
	-cd Cufflinks/ -cc mm9-EnsGene -co CoverageFiles/mm9-EnsGene.coverage \
	-ch chr1 -o Optimizer-EnsGene -t Cufflinks/mm9-EnsGene.transcripts.gtf.tmap

The results of this analysis can be found in the directory ./Optimizer-EnsGene-chr1/ as the Optimizer script concatenates the chromosome used to the specified output directory. TSS-Optimizer produces a number of files including a chromosome.coverage file, which is simply a reduced version of the .coverage file. TSS-Optimizer.pl produces all of the files produced by TSS-Predictor.pl (see section 4). Multiple copies of isotigs.trimmed.gtf, distances.txt and tss.prediction.bed are produced, one for each cutoff fraction interrogated. 

The most informative output of TSS-Optimizer.pl is the graph Optimizer-MADPlot.jpg (Figure 1). This barchart shows the MAD between predicted TSSs and the nearest annotated TSS for each gene over all cutoff fractions assessed. The cutoff fraction that gives the lowest MAD is shown in red, this should be given to TSS-Predictor.pl to identify TSS usage genome wide. The actual MAD calculated by TSS-Optimizer contains no post-filtering (see section 3.7) and is typically much greater than the value obtained in the final predictions.

![Barchart showing the median absolute deviation between predicted transcription start sites and the nearest annotated TSS for the gene over the coverage cutoff fractions assessed by TSS-Optimizer.pl. The cutoff fraction with the minimum MAD is shown in red (0.04).](./Optimizer-MADPlot.pdf)



3.6. Running TSS-Predictor.pl
-----------------------------

In this next step we run TSS-Predictor to predict transcription start sites genome wide. This script attempts to identify the major transcription start site for each reference gene assigned to a a Cufflinks gene_id by cuffcompare. The output will include erroneous TSS calls due to a number of factors, these include incorrect transcript assembly and selection of the wrong Cufflinks transfrag for trimming. Strategies exist to get high confidence TSS predictions, see next section. 
	
	
	TSS-Predictor.pl -a <Annotations.gtf> -cd <Cufflinks_directory> \
	-cc <Cuffcompare_prefix> -t <tmap_file> -co <Coverage_file>

+---------+---------------------+--------------------------------------------+
| Options | Default             | Description                                |
+=========+=====================+============================================+
| -a      | []                  | The main reference annotation file         |
|         |                     | i.e. the one used with Cuffcompare         |
+---------+---------------------+--------------------------------------------+
| -cd     | []                  | Path to the Cufflinks directory            |
+---------+---------------------+--------------------------------------------+
| -cc     | []                  | Cuffcompare prefix                         |
+---------+---------------------+--------------------------------------------+
| -t      | []                  | Path to Cuffcompare transcripts.gtf.tmap   |
+---------+---------------------+--------------------------------------------+
| -co     | []                  | Path to the .coverage file produced by     |
|         |                     | GetCoverage.pl                             |
+---------+---------------------+--------------------------------------------+
| -f      | [0.1]               | Cutoff fraction (of expected per-base      |
|         |                     | coverage) to trim annotations              |
+---------+---------------------+--------------------------------------------+
| -u      | [1]                 | Number of bases upstream of predicted TSS  |
|         |                     | to report                                  |
+---------+---------------------+--------------------------------------------+
| -d      | [0]                 | Number of bases downstream of predicted    |
|         |                     | TSS to report                              |
+---------+---------------------+--------------------------------------------+


Next, predict TSS usage genome-wide with a coverage cutoff of 0.04:-

	cd ~/TSS-Predictor-Example
	
	TSS-Predictor.pl -f 0.04 -a Annotations/Mus_musculus.mm9.67.gtf \
	-cd Cufflinks/ -cc mm9-EnsGene -co CoverageFiles/mm9-EnsGene.coverage \
	-t Cufflinks/mm9-EnsGene.transcripts.gtf.tmap



3.7. Post-filtering for sensible results
----------------------------------------

TSS-Predictor.pl produces a few plots to help you choose an appropriate method for cleaning up your data. The final step in the TSS-Predictor pipeline is SelectTSS.pl. This takes user-defined cutoffs and reports only transcription start sites that pass these filter(s). For SelectTSS.pl usage see the next section.

###Expression-based filtering####
Transcripts that are expressed at a low level are less likely to be assembled correctly, and hence TSS predictions for the least abundant transcripts are less likely to be accurate. Expression-based filtering involves disregarding genes expressed below a user-defined level, optimal cutoffs vary between experiments and depend on a number of factors including sequencing depth and read length. In some cases expression-based filtering may not be suitable at all. TSS-Predictor uses R to produce two plots to aid the selection of appropriate cutoffs. 

Figure 2 shows the MAD expression plot produced by TSS-Optimizer for our sample dataset, this clearly demonstrates the effect of gene expression on TSS prediction accuracy. We assume that one of the annotated transcription start sites for each expressed gene is correct. For high quality, well annotated genomes such as mouse and human this appears to be a relatively safe assumption. The MAD expression plot shows the median absolute deviation for genes, grouped into bins based on the per-base coverage (expression). We observe a much higher variation in the distances between predicted and annotated TSS for genes expressed at a lower level. 

![The MAD expression plot.](./MadExpressionPlot.pdf)

The dashed black line in Figure 2 shows the MAD for all genes in the dataset.  TSS-Predictor advises an expression cutoff as the median per-base gene expression of the first bin below the line. The effect of trimming at this cutoff can be seen in the ExprCutoffPlot.jpg file produced by TSS-Predictor.pl, which for our example dataset is shown in Figure 3. This shows that a large number of lowly expressed outliers will be removed by disreagarding genes expressed below the red line. In principle this approach facilitates the discovery of novel TSSs by identifying highly expressed genes with predicted TSS a large distance from the nearest predicted TSS for the gene. In practise however, these 'novel' TSSs most often arise from incorrect transcript assembly, or selection of an inappropriate reference transfrag to trim. For some datasets expression-based filtering does not work well, usually because the distinction between low expression (and high variance) and high expression (and low variance) is not particuarly clear. In this case an alternative approach to filtering is distance based, see below. 

![The expression cutoff plot. TSS-Predictor reports the value of the expression cutoff to the command line. In this instance it is 65.09, this value is printed to the terminal during TSS-Predictor analysis.](./ExprCutoffPlot.pdf)


###Distance based filtering###
Distance based filtering assumes explicitly that one of the annotated transcription start sites for each gene is correct. Genes with a TSS greater than a given distance from the nearest predicted TSS are disregarded. 500 bases often appears to be a good cutoff. The DistanceCutoffPlot (Figure 4) shows the effect of removing TSS greater than 500bp from the nearest annotated TSS for the gene.

![The distance cutoff plot.](./DistanceCutoffPlot.pdf)


3.8. Selecting a high confidence TSS set using SelectTSS.pl
-----------------------------------------------------------

SelectTSS.pl is used to filter TSS predictions. A tss.selected.bed file is produced which contains only transcription start sites that pass either expression-based or distance-based filtering. If both a distance and expression cutoff are supplied then only transcripts that pass both filters are reported. There is an additional option to specify whether to report the predicted TSS, or the TSS of the nearest annotated transcript for the gene.

	SelectTSS.pl [options -d -e -r] -o <output prefix>

+---------+---------------------+--------------------------------------------+
| Options | Default             | Description                                |
+=========+=====================+============================================+
| -o      | tss.selected        | Output file prefix                         |
+---------+---------------------+--------------------------------------------+
| -d      | []                  | TSS distance cutoff                        |
+---------+---------------------+--------------------------------------------+
| -e      | []                  | Expression level cutoff                    |
+---------+---------------------+--------------------------------------------+
| -r      | []                  | Report nearest reference transcription     |
|         |                     | start site                                 |
+---------+---------------------+--------------------------------------------+

We want to select only high confidence TSS's in our example. The following command will report predicted TSS's that satisfy both distance AND expression criteria:- 

	cd ~/...
	SelectTSS.pl -d 500 -e 65.09 


4.0. Output files and formats
=============================

4.1. Coverage file
------------------

The coverage file produced by GetCoverage.pl is a tab-delimited text file containing the per-base coverage of each transcribed location on the reference genome as determined by Cufflinks. The file format is as follows:-

+---------------+---------------+----------------------------------------------+
|Field position |Field name     |Desription                                    |
+===============+===============+==============================================+
|1              |Chromosome     |Reference chromosome/scaffold                 |
+---------------+---------------+----------------------------------------------+
|2              |Position       |Position on the reference chromosome          |
+---------------+---------------+----------------------------------------------+
|3              |Read depth     |Read coverage at position on genome           |
+---------------+---------------+----------------------------------------------+


4.2. trim.combined.gtf
----------------------

The trim.combined.gtf file is produced by Prep4Trim.pl. This GTF file is a reduced version of the Cuffcompare combined.gtf file and contains one transfrag record for each expressed reference gene (as determined by Cuffcompare). Specifically, the longest transfrag on the same strand as the reference gene is reported. For a detailed descripion combined.gtf format see the Cufflinks manual <http://cufflinks.cbcb.umd.edu/manual.html>.

4.3. trim.gtf.tmap
------------------

The trim.gtf.tmap file is an amended version of the transcripts.gtf.tmap file produced by Cuffcompare. This tab seperated file contains only the reference isotigs selected for trimming. There are also a couple of changes in this file relative to the Cuffcompare output. Firstly, the _ref\_id_ field contains the reference gene, as defined by the _gene\_id_ tag in the annotations GTF file, rather than the reference transcript. The second difference is the _gene.cov_ field. The cuffcompare output contains a _cov_ field reporting the expected transcript coverage, whereas the trim.gtf.tmap contains the expected per-base coverage for the gene calculated by summing the length-normalized coverage estimates for all transcripts in the gene (see section 5.1). 

TSS-Predictor makes use of fields 2 and 10, the rest of the file is produced by Cufflinks. The fields in the trim.gtf tmap file are as follows:

+---------------+---------------+----------------------------------------------+
|Field position |Field name     |Desription                                    |
+===============+===============+==============================================+
|1              |ref\_gene\_id  |Gene name assigned by Cuffcompare             |
|               |               |(_gene\_name_ attribute of GTF file)          |
+---------------+---------------+----------------------------------------------+
|2              |ref\_id        |_gene\_id_ attribute of the the GTF file      |
+---------------+---------------+----------------------------------------------+
|3              |class\_code    | Cuffcompare class code                       |
+---------------+---------------+----------------------------------------------+
|4              |cuff\_gene\_id | The gene id assigned by Cufflinks. Multiple  |
|               |               |_cuff\_gene\_id_'s may map to a single        |
|               |               |_ref\_gene\_id_ attribute                     |
+---------------+---------------+----------------------------------------------+
|5              |cuff\_id       | Name of the Cufflinks transcript id selected |
|               |               | for trimming for _ref\_gene\_id_             |
+---------------+---------------+----------------------------------------------+
|6              |FMI            | Expression of this isoform relative to major | 
|               |               | isoform for the gene (defined by Cufflnks,   |
|               |               | not used)                                    |
+---------------+---------------+----------------------------------------------+
|7              |FPKM           | FPKM for the transcript _cuff\_id_. Note this|
|               |               | is the transcript (NOT gene) expression      |
+---------------+---------------+----------------------------------------------+
|8              |FPKM\_conf\_lo | Upper confidence interval for FPKM           |
+---------------+---------------+----------------------------------------------+
|9              |FPKM\_conf\_hi | Lower confidence interval for FPKM           |
+---------------+---------------+----------------------------------------------+
|10             |gene.cov       | The per-base coverage estimate for the gene  |
|               |               | as used by TSS-Predictor                     |
+---------------+---------------+----------------------------------------------+
|11             |len            | _cuff\_id_ transcript length                 |
+---------------+---------------+----------------------------------------------+
|12             |major\_iso\_id | _cuff\_id_ of the genes (_cuff\_gene\_id_)   |
|               |               | major isoform                                |
+---------------+---------------+----------------------------------------------+


4.4. isotigs.trimmed.gtf
------------------------

The isotigs.trimmed.gtf file is an amended version of the trim.combined.gf file (see section 4.1). An additional flag, _contain\_tss_, is added to field 9, this contains either "yes" or "no" and identifies the exons which are predicted to contain the primary TSS for the gene. For transcripts on the '+' strand with _contain\_tss_ = "yes" the exon start position (field 4) is amended to reflect the predicted TSS. When _contain\_tss_ = "yes" and the transcript is on the '-' strand the exon end position (field 5) is amended accordingly.


4.5. tss.prediction.bed
-----------------------

This BED (tab delimited) file contains the TSS predictions for each expressed reference gene. By default the TSS lenth is 1 base. Proximal promoter regions flanking the predicted TSS can be extracted by providing TSS-Predictor with the _-u_ and _-d_ flags. Please note that changing the values of _-u_ and _-d_ will affect the distance estimates produced in the distances.txt file. The format of tss.predictions.bed is as follows:-

+---------------+---------------+---------------------------------+
|Field position |Field name     |Desription                       |
+---------------+---------------+---------------------------------+
|1              |Chromosome     | Reference chromosome/scaffold   |
+---------------+---------------+---------------------------------+
|2              |Start          | TSS prediction start            |
+---------------+---------------+---------------------------------+
|3              |End            | TSS prediction end              | 
+---------------+---------------+---------------------------------+
|4              |Gene name      | The reference gene id           | 
+---------------+---------------+---------------------------------+
|5              |Score          | Not used                        |
+---------------+---------------+---------------------------------+
|6              |Strand         | Strand the reference gene is on |
+---------------+---------------+---------------------------------+


4.6. distances.txt
------------------

The tab-delimited distances.txt file contains a variety of information for each expressed reference gene. The format of this file is as follows:-

+---------------+---------------------------+----------------------------------+
|Field position |Field name   	            |Desription                        |
+---------------+---------------------------+----------------------------------+
|1              |Gene                       |Reference gene id                 |
+---------------+---------------------------+----------------------------------+
|2              |Chr                        |Chromosome                        |
+---------------+---------------------------+----------------------------------+
|3              |Prediction\_start          |Predicted TSS                     |
+---------------+---------------------------+----------------------------------+
|4              |Nearest\_transcript        |Closese reference transcript to   |
|               |                           |TSS prediction                    |
+---------------+---------------------------+----------------------------------+
|5              |Nearest\_transcript\_start |TSS of nearest annotated reference|
|               |                           |transcript                        |
+---------------+---------------------------+----------------------------------+
|6              |Distance                   |Distance between predicted and    |
|               |                           |nearest reference TSS             |
+---------------+---------------------------+----------------------------------+


4.7. tss.selected.bed
---------------------

This BED file contains the subset of tss.prediction.bed containing high-confidence TSS predictions (i.e. those which have passed the filters set by SelectTSS.pl). 


-------------------

5.0. Overview of the algorithm 
==============================

TSS-Predictor works by selecting and trimming Cufflinks gene annotations based upon the ratio between the genes expected coverage (expression level) and per-base coverage at each position on the gene. Below is an overview of the steps involved:-

1. Firstly, reads are mapped to the reference genome allowing multi-mapping reads. Because the Cufflinks output is used by TSS-Predictor I use TopHat, but this is not essential. **Depth is important for TSS prediction, for this reason I merge BAM alignments after read mapping**. 

2. Cufflinks is then run with the default parameters and transcribed fragments 'transfrags' are linked to reference genes using Cuffcompare. Ensembl is my favourite annotation source, but all you need is a many-to-one transcript\_id to gene\_id relationship in the GTF file. Note UCSC GTF files obtained from the table browser do not maintain this relationship.

3. A coverage file is produced using _GetCoverage.pl_. This script makes use of the Cufflinks output and BEDTools to calculate the per-base read coverage for every transcribed base reported in the Cufflinks cuffcompare.combined.gtf file. 

4. TSS-Optimizer uses a reduced dataset of one chromosome to determine optimal conditions (see Section 3.5 _Using TSS-Optimizer to identify optimal coverage cutoffs_). TSS-Predictor predicts major TSS utilization genome-wide once the optimal conditions have been determined. TSS-Optimizer and TSS-Predictor run the same core algorithm and are described together below:-

5.1. Prep4Trim.pl
-----------------

Prep4Trim.pl calculates the expected coverage for each reference gene and selects one representative transfrag for trimming. Expected gene expression is calculated by summing length-normalized expression (coverage) for all Cufflinks transfrags associated with each reference gene, as defined by the _gene\_id_ tag in the reference GTF file (Figure 5). Additionally, one representative transfrag for each reference gene is selected for trimming. Specifically, the longest Cufflinks transfrag on the same strand as the reference gene is selected.  

![Gene-wise coverage calculation. Each expressed reference gene _G_ is composed of _n_ transfrags (_T_), each with coverage (_C_) and length (_L_). Thus, _C_~_T__i_~  is the coverage for the _i_^th^ transfrag of gene _G_, and _L_~_T__i_~ is the length of the _i_^th^ transfrag of gene _G_.](./Normalization.pdf )

Prep4Trim.pl reports two files - trim.combined.gtf and trim.gtf.tmap. Trim.combined.gtf is a reduced version of the Cufflinks cuffcompare.combined.gtf file, and contains only the transfrags selected for trimming in the next stage of the algorithm. Trim.gtf.tmap is an amended version of the cuffcompare .tmap output. This contains a mapping between each Cufflinks transfrag to be trimmed and the reference gene, as well as the amended coverage estimate for the gene (see section 4.3). All other fields in trim.gtf.tmap file are transfrag level information generated by Cuffcompare.


5.2. Trim\_Cufflinks\_GTF.pl
----------------------------

The core of the algorithm. Trim\_Cufflinks\_GTF.pl takes the expected gene-level coverage information from trim.gtf.tmap and transfrag information from trim.combined.gtf. Using the per-base coverage information calculated by _GetCoverage.pl_ each selected isotig is trimmed at a fraction of the expected coverage for the gene (_X_). 

The algorithm iterates, one base at a time, from the last base (3') in the first exon of each transfrag (providing it is covered by >_X_ * _C_~_G_~) towards the first base. The GTF file is trimmed when the per-base coverage of transfrag _T_~_max_~ drops below _X_ * _C_~_G_~. If the last base of exon 1 is expressed at a lower level than _X_ then the second exon is examined, so on and so forth...

The value of _X_ is provided by the user, and determined through empirical optimization. This is what the TSS-Optimizer is for. Optimal values for _X_ are then used to predict transcription start sites genome-wide using TSS-Predictor.pl.

Trim\_Cufflinks\_GTF.pl reports an isotigs.trimmed.gtf file. This is an amended version of trim.combined.gtf with revised transcription start site information. An additional variable _"contain\_tss"_ is added to the 9^th^ field of the GTF file to identify exons predicted to contain transcription start sites.


5.3. Get-TSS.pl
---------------

This script links reports a predicted transcription start site for each reference gene in BED format. By default this is tss.prediction.bed. The _-u_ and _-d_ flags are used to set the number of bases upstream (_u_) and downstream (_d_) of the predicted start site to report.

5.4. GetDistances.pl
--------------------

GetDistances.pl generates the distances.txt file. This tab-delimited file contains the reference gene, chromosome, predicted transcription start site, nearest reference transcript, nearest reference transcription start site and the distance between the two predictions.

5.5. SelectTSS.pl
-----------------

SelectTSS.pl acts as a filter. It reporting only isoforms which are expressed above an expression cutoff defined by _-e_ or a distance cutoff defined by _-d_. The distance metric used is the distance between the predicted TSS and the nearest annotated TSS for the associated gene. SelectTSS.pl makes use of the trim.gtf.tmap and distances.txt files to obtain expression and distance information respectively and reports two files tss.selected.bed amd tss.selected.distances.txt. These files are identical in format to tss.prediction.bed and distances.txt.

5.6. CollateDistances.pl (Optimizer only)
-----------------------------------------

This script is run called by TSS-Optimizer and produces the CutoffDistances.txt file used to determine the optimal trimming conditions.

