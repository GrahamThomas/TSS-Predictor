<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
  <meta http-equiv="Content-Style-Type" content="text/css" />
  <meta name="generator" content="pandoc" />
  <title></title>
</head>
<body>
<h1 id="tss-predictor">TSS-Predictor</h1>
<p>Version 1.0 by Graham Thomas</p>
<hr />
<h1 id="contents">Contents</h1>
<h3 id="introduction">1.0. Introduction</h3>
<h3 id="getting-started">2.0. Getting Started</h3>
<pre><code>2.1. Dependencies
2.2. Installation</code></pre>
<h3 id="using-tss-predictor">3.0. Using TSS-Predictor</h3>
<pre><code>3.1. Preparing your data
3.2. Mapping reads to the reference genome
3.3. Running Cufflinks and Cuffcompare
3.4. Generating a .coverage file with GetCoverage.pl
3.5. Using TSS-Optimizer to identify optimal coverage cutoffs
3.6. Running TSS-Predictor.pl
3.7. Post-filtering for sensible results
3.8. Seecting a high confidence TSS set using SelectTSS.pl</code></pre>
<h3 id="output-files-and-formats">4.0. Output files and formats</h3>
<pre><code>4.1. Coverage file
4.2. trim.combined.gtf
4.3. trim.gtf.tmap
4.4. isotigs.trimmed.gtf
4.5. tss.prediction.bed
4.6. distances.txt
4.7. tss.selected.bed</code></pre>
<h3 id="overview-of-the-algorithm">5.0. Overview of the algorithm</h3>
<pre><code>5.1. Prep4Trim.pl
5.2. Trim_Cufflinks_GTF.pl
5.3. Get-TSS.pl
5.4. GetDistances.pl
5.5. SelectTSS.pl
5.6. CollateDistances.pl</code></pre>
<h3 id="references">6.0. References</h3>
<hr />
<h1 id="introduction-1">1.0 Introduction</h1>
<p>TSS-Predictor has been developed to identify promoter usage for genes in RNA-Seq experiments. It is designed to be used when both a high quality reference genome and annotation set are available, for example human or mouse. The output of TSS-Predictor is a BED file containing predicted transcription start sites for each reference gene, these may then be used for downstream applications such as <em>cis</em>-regulatory analysis.</p>
<p>TSS-Predictor is written in perl and is available at GitHub <a href="https://github.com/GrahamThomas/TSS-Predictor"><code class="url">https://github.com/GrahamThomas/TSS-Predictor</code></a> and as an online tool at Galaxy (Not Yet!!!) and GeneProf (Not Yet!!!).</p>
<h1 id="getting-started-1">2.0. Getting Started</h1>
<h2 id="dependencies">2.1. Dependencies</h2>
<p>In order to use TSS-Predictor you must have working versions of the following in your $PATH:-</p>
<ol style="list-style-type: decimal">
<li>A short read mapper compatible with Cufflinks. I use TopHat <a href="http://tophat.cbcb.umd.edu/"><code class="url">http://tophat.cbcb.umd.edu/</code></a>.</li>
<li>Cufflinks <a href="http://cufflinks.cbcb.umd.edu/"><code class="url">http://cufflinks.cbcb.umd.edu/</code></a></li>
<li>BEDTools <a href="http://code.google.com/p/bedtools/"><code class="url">http://code.google.com/p/bedtools/</code></a> and</li>
<li>R, with the libraries ggplot2 and reshape2 installed <a href="http://cran.r-project.org/"><code class="url">http://cran.r-project.org/</code></a></li>
</ol>
<p>To follow this guide you will need to merge and sort some BAM files. For this you will also need SAMtools <a href="http://samtools.sourceforge.net/"><code class="url">http://samtools.sourceforge.net/</code></a>.</p>
<h2 id="installation">2.2. Installation</h2>
<p>TSS-Predictor is implemented in Perl and R. To get up and running all you need to do is:-</p>
<ol style="list-style-type: decimal">
<li><p>Add TSS-Predictor/bin to your $PATH variable:-</p>
<p>$PATH=$PATH:/&lt;Path to installation directory&gt;/TSS-Predictor/bin</p></li>
<li><p>And give executable permission to All scripts in the directory.</p>
<p>cd &lt;Path to installation directory&gt;/TSS-Predictor/bin</p>
<p>chmod +x *</p></li>
</ol>
<hr />
<h1 id="using-tss-predictor-1">3.0. Using TSS-Predictor</h1>
<h2 id="preparing-your-data">3.1. Preparing your data</h2>
<p>Below is a step-by-step guide on how to use TSS-Predictor beginning with raw illumina reads. In this example we use 51 base paired-end Illumina RNA-Seq data from a macrophage transcriptomics project (REF). Raw data and BAM files for this project are available from the SRA (ACCESSION... <a href="http://www.ebi.ac.uk/ena/"><code class="url">http://www.ebi.ac.uk/ena/</code></a>). All of the intermediate files generated in the working example are available at our FTP site <link...>.</p>
<pre><code>cd ~/TSS-Predictor-Example/RawData/
ls
4Ne_1_F.fastq  4Ne_3_R.fastq  4TG_3_F.fastq  BNe_2_R.fastq  BTG_2_F.fastq            
4Ne_1_R.fastq  4TG_1_F.fastq  4TG_3_R.fastq  BNe_3_F.fastq  BTG_2_R.fastq  
4Ne_2_F.fastq  4TG_1_R.fastq  BNe_1_F.fastq  BNe_3_R.fastq  BTG_3_F.fastq  
4Ne_2_R.fastq  4TG_2_F.fastq  BNe_1_R.fastq  BTG_1_F.fastq  BTG_3_R.fastq
4Ne_3_F.fastq  4TG_2_R.fastq  BNe_2_F.fastq  BTG_1_R.fastq</code></pre>
<p>A pre-built Bowtie index for mouse (mm9) was downloaded from <a href="http://bowtie-bio.sourceforge.net/index.shtml"><code class="url">http://bowtie-bio.sourceforge.net/index.shtml</code></a>.</p>
<pre><code>cd ~/TSS-Predictor-Example/BowtieIndexes/
wget ftp://ftp.cbcb.umd.edu/pub/data/bowtie_indexes/mm9.ebwt.zip
unzip mm9.zip</code></pre>
<p>Reconstitute a reference genome from the bowtie index.</p>
<pre><code>cd ~/TSS-Predictor-Example/BowtieIndexes/
bowtie-inspect mm9 &gt; mm9.fa</code></pre>
<p>An Ensembl GTF file containing reference annotations was downloaded from the Ensembl FTP site <a href="ftp://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/"><code class="url">ftp://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/</code></a>. Ensembl annotations are preferable as, unlike UCSC annotations, a many-to-one relationship exists between gene 'gene_id' and 'transcript_id' fields of the GTF file. This relationship is required to combine per-transcript expression levels to per-gene expression levels.</p>
<pre><code>cd ~/TSS-Predictor-Example/Annotations/ 
wget ftp://ftp.ensembl.org/pub/release-67/gtf/\
mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz </code></pre>
<p>Ensembl chromosomes need converting from NCBIM37 to mm9 style. This can be done with Ensembl2UCSC.pl. This script prints the annotations in the same order as the FASTA headers in the supplied genome and removes annotations belonging to the 'random' chromosomes. An alternative approach would be to map reads directly to the Ensembl genome build, thereby ommitting this step.</p>
<pre><code>Ensembl2UCSC.pl -r ../BowtieIndexes/mm9.fa -e Mus_musculus.NCBIM37.67.gtf \ 
-o Mus_musculus.mm9.67.gtf</code></pre>
<h2 id="mapping-reads-to-the-reference-genome">3.2. Mapping reads to the reference genome</h2>
<p>Map reads to the reference genome using TopHat (version v2.0.0, Bowtie version 0.1.18.0 were used in this example). Multi-mapping reads are allowed by using the default parameters. Multi-mapped reads are important as they facilitate the assembly of transcripts containing repetitive sequence. In downstream analysis it may be wise to disregard transcripts containing a large proportion of multi-mapped reads as they may be unreliable.</p>
<p>Transcription start site prediction is more effective with greater sequencing depth. My preferred strategy is to map reads for each lane individually, and then merge the alignments prior to transcript assembly with Cufflinks.</p>
<p>Map the reads with TopHat:-</p>
<pre><code>cd ~/TSS-Predictor-Example/TopHat
tophat -p 5 -o KNe_1_multi_out ../BowtieIndexes/mm9 \
../RawData/4Ne_1_F.fastq ../RawData/4Ne_1_R.fastq
...
tophat -p 5 -o BTG_3_multi_out ../BowtieIndexes/mm9 \
../RawData/BTG_3_F.fastq ../RawData/BTG_3_R.fastq</code></pre>
<p>Merge and sort the alignments. I use SAMtools <a href="http://samtools.sourceforge.net/"><code class="url">http://samtools.sourceforge.net/</code></a>:-</p>
<pre><code>cd ~/TSS-Predictor-Example/TopHat/

###Reconstitute a SAM header from one of the TopHat alignments
samtools view -H BNe_1_multi_out/accepted_hits.bam &gt; mm9.header.sam


samtools cat -h mm9.header.sam -o Macrophage.cat.bam \
BNe_1_multi_out/accepted_hits.bam BNe_2_multi_out/accepted_hits.bam \
BNe_3_multi_out/accepted_hits.bam KNe_1_multi_out/accepted_hits.bam \ 
KNe_2_multi_out/accepted_hits.bam KNe_3_multi_out/accepted_hits.bam \
BTG_1_multi_out/accepted_hits.bam BTG_2_multi_out/accepted_hits.bam \
BTG_3_multi_out/accepted_hits.bam KTG_1_multi_out/accepted_hits.bam \
KTG_2_multi_out/accepted_hits.bam KTG_3_multi_out/accepted_hits.bam

samtools sort Macrophage.cat.bam Macrophage.sorted</code></pre>
<h2 id="running-cufflinks-and-cuffcompare">3.3. Running Cufflinks and Cuffcompare</h2>
<p>Assemble the mapped reads into transcribed fragments 'transfrags' using Cufflinks (here we have used version 2.0.2). In this example we have not provided reference annotations to guide transcript assembly. You may do so if you wish.</p>
<pre><code>cd ~/TSS-Predictor-Example
cufflinks -p 8 -L Macrophage-Cufflinks -o Cufflinks TopHat/Macrophage.sorted.bam</code></pre>
<p>Assign the transfrags to reference gene models using Cuffcompare</p>
<pre><code>cd ~/TSS-Predictor-Example/Cufflinks
cuffcompare -r ../Annotations/Mus_musculus.mm9.67.sorted.gtf -M \
-s ../mm9.fa -o mm9-EnsGene transcripts.gtf</code></pre>
<h2 id="generating-a-.coverage-file-with-getcoverage.pl">3.4. Generating a .coverage file with GetCoverage.pl</h2>
<p>GetCoverage.pl makes use of BEDTools to calculate the per-base read coverage at the genomic locations reported in the combined.gtf file produced by Cuffcompare. The script takes the following arguments:-</p>
<pre><code>GetCoverage.pl -c &lt;cuffcompare.combined.gtf&gt; -b &lt;mapped_reads.bam&gt; \
-g &lt;.genome file&gt; -o &lt;output&gt;</code></pre>
<table>
<col width="13%" />
<col width="62%" />
<thead>
<tr class="header">
<th align="left">Options</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">-c</td>
<td align="left">Path to the cuffcompare combined.gtf file</td>
</tr>
<tr class="even">
<td align="left">-b</td>
<td align="left">Path to the mapped reads BAM file</td>
</tr>
<tr class="odd">
<td align="left">-g</td>
<td align="left">Path to the .genome file</td>
</tr>
<tr class="even">
<td align="left">-o</td>
<td align="left">Name of the output .coverage file</td>
</tr>
</tbody>
</table>
<p>The .genome file is a two column tab-delimited file stating each chromosome/scaffold name and it's size in base pairs. BEDTools is shipped with .genome files for mouse and human reference genomes. If you are working with a different species you will have to create one yourself, see the BEDTools manual (<a href="http://code.google.com/p/bedtools/"><code class="url">http://code.google.com/p/bedtools/</code></a>) for more information.</p>
<p>Making the .coverage file can take some time. We make the .coverage file for this example as follows:-</p>
<pre><code>cd ~/TSS-Predictor-Example/CoverageFiles
cp &lt;Path To BEDTools&gt;/genomes/mouse.mm9.genome .

GetCoverage.pl -c ../Cufflinks/mm9-EnsGene.combined.gtf \
-b ../TopHat/Macrophage.sorted.bam -g mouse.mm9.genome \
-o mm9-EnsGene.coverage</code></pre>
<h2 id="using-tss-optimizer.pl-to-identify-optimal-cutoffs">3.5. Using TSS-Optimizer.pl to identify optimal cutoffs</h2>
<p>TSS-Optimizer.pl is a wrapper script for the main TSS prediction algorithm. For a description of how the main algorithm works see section 5 or the supplementary information of (REF MY PAPER). TSS-Optimizer.pl takes one chromosome, defined by the user and predicts transcription start sites over a range of cutoff fractions. TSS-Optimizer.pl reports, for each gene at each cutoff fraction, the distance between the predicted TSS and nearest annotated TSS for the same gene. TSS-Optimizer also produces a graph showing the median absolute deviation (MAD) between the predicted and nearest annotated transcription start sites at each cutoff fraction. We choose the cutoff fraction that gives the lowest MAD for use with TSS-Predictor.</p>
<pre><code>TSS-Optimizer.pl -a &lt;Annotations.gtf&gt; -cd &lt;Cufflinks_directory&gt; \
-cc &lt;Cuffcompare_prefix&gt; -t &lt;tmap_file&gt; -co &lt;Coverage_file&gt; -ch &lt;chromosome&gt;</code></pre>
<table>
<col width="12%" />
<col width="20%" />
<col width="66%" />
<thead>
<tr class="header">
<th align="left">Options</th>
<th align="left">Default</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">-a</td>
<td align="left">[]</td>
<td align="left">The main reference annotation file i.e. the one used with Cuffcompare</td>
</tr>
<tr class="even">
<td align="left">-cd</td>
<td align="left">[]</td>
<td align="left">Path to the Cufflinks directory</td>
</tr>
<tr class="odd">
<td align="left">-cc</td>
<td align="left">[]</td>
<td align="left">Prefix used when running Cuffcompare</td>
</tr>
<tr class="even">
<td align="left">-t</td>
<td align="left">[]</td>
<td align="left">Path to the Cuffcompare transcripts.gtf.tmap file</td>
</tr>
<tr class="odd">
<td align="left">-co</td>
<td align="left">[]</td>
<td align="left">Path to the .coverage file produced by GetCoverage.pl</td>
</tr>
<tr class="even">
<td align="left">-ch</td>
<td align="left">[]</td>
<td align="left">Name of chromosome for optimization</td>
</tr>
<tr class="odd">
<td align="left">-s</td>
<td align="left">[0.01]</td>
<td align="left">Coverage cutoff fraction to begin trimming</td>
</tr>
<tr class="even">
<td align="left">-e</td>
<td align="left">[0.20]</td>
<td align="left">Coverage cutoff fraction to end trimming</td>
</tr>
<tr class="odd">
<td align="left">-i</td>
<td align="left">[0.01]</td>
<td align="left">Increment, or step size, between lower and upper coverage cutoff fractions</td>
</tr>
<tr class="even">
<td align="left">-o</td>
<td align="left">[./Optimizer]</td>
<td align="left">Output directory (this will be suffixed with the selected chromosome)</td>
</tr>
</tbody>
</table>
<p>Running TSS-Optimizer for our data set is simple. We use the following command:-</p>
<pre><code>cd ~/TSS-Predictor-Example

TSS-Optimizer.pl -a Annotations/Mus_musculus.mm9.67.gtf \
-cd Cufflinks/ -cc mm9-EnsGene -co CoverageFiles/mm9-EnsGene.coverage \
-ch chr1 -o Optimizer-EnsGene -t Cufflinks/mm9-EnsGene.transcripts.gtf.tmap</code></pre>
<p>The results of this analysis can be found in the directory ./Optimizer-EnsGene-chr1/ as the Optimizer script concatenates the chromosome used to the specified output directory. TSS-Optimizer produces a number of files including a chromosome.coverage file, which is simply a reduced version of the .coverage file. TSS-Optimizer.pl produces all of the files produced by TSS-Predictor.pl (see section 4). Multiple copies of isotigs.trimmed.gtf, distances.txt and tss.prediction.bed are produced, one for each cutoff fraction interrogated.</p>
<p>The most informative output of TSS-Optimizer.pl is the graph Optimizer-MADPlot.jpg (Figure 1). This barchart shows the MAD between predicted TSSs and the nearest annotated TSS for each gene over all cutoff fractions assessed. The cutoff fraction that gives the lowest MAD is shown in red, this should be given to TSS-Predictor.pl to identify TSS usage genome wide. The actual MAD calculated by TSS-Optimizer contains no post-filtering (see section 3.7) and is typically much greater than the value obtained in the final predictions.</p>
<div class="figure">
<embed src="./Optimizer-MADPlot.pdf"><p class="caption">Barchart showing the median absolute deviation between predicted transcription start sites and the nearest annotated TSS for the gene over the coverage cutoff fractions assessed by TSS-Optimizer.pl. The cutoff fraction with the minimum MAD is shown in red (0.04).</p>
</div>
<h2 id="running-tss-predictor.pl">3.6. Running TSS-Predictor.pl</h2>
<p>In this next step we run TSS-Predictor to predict transcription start sites genome wide. This script attempts to identify the major transcription start site for each reference gene assigned to a a Cufflinks gene_id by cuffcompare. The output will include erroneous TSS calls due to a number of factors, these include incorrect transcript assembly and selection of the wrong Cufflinks transfrag for trimming. Strategies exist to get high confidence TSS predictions, see next section.</p>
<pre><code>TSS-Predictor.pl -a &lt;Annotations.gtf&gt; -cd &lt;Cufflinks_directory&gt; \
-cc &lt;Cuffcompare_prefix&gt; -t &lt;tmap_file&gt; -co &lt;Coverage_file&gt;</code></pre>
<table>
<col width="12%" />
<col width="28%" />
<col width="58%" />
<thead>
<tr class="header">
<th align="left">Options</th>
<th align="left">Default</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">-a</td>
<td align="left">[]</td>
<td align="left">The main reference annotation file i.e. the one used with Cuffcompare</td>
</tr>
<tr class="even">
<td align="left">-cd</td>
<td align="left">[]</td>
<td align="left">Path to the Cufflinks directory</td>
</tr>
<tr class="odd">
<td align="left">-cc</td>
<td align="left">[]</td>
<td align="left">Cuffcompare prefix</td>
</tr>
<tr class="even">
<td align="left">-t</td>
<td align="left">[]</td>
<td align="left">Path to Cuffcompare transcripts.gtf.tmap</td>
</tr>
<tr class="odd">
<td align="left">-co</td>
<td align="left">[]</td>
<td align="left">Path to the .coverage file produced by GetCoverage.pl</td>
</tr>
<tr class="even">
<td align="left">-f</td>
<td align="left">[0.1]</td>
<td align="left">Cutoff fraction (of expected per-base coverage) to trim annotations</td>
</tr>
<tr class="odd">
<td align="left">-u</td>
<td align="left">[1]</td>
<td align="left">Number of bases upstream of predicted TSS to report</td>
</tr>
<tr class="even">
<td align="left">-d</td>
<td align="left">[0]</td>
<td align="left">Number of bases downstream of predicted TSS to report</td>
</tr>
</tbody>
</table>
<p>Next, predict TSS usage genome-wide with a coverage cutoff of 0.04:-</p>
<pre><code>cd ~/TSS-Predictor-Example

TSS-Predictor.pl -f 0.04 -a Annotations/Mus_musculus.mm9.67.gtf \
-cd Cufflinks/ -cc mm9-EnsGene -co CoverageFiles/mm9-EnsGene.coverage \
-t Cufflinks/mm9-EnsGene.transcripts.gtf.tmap</code></pre>
<h2 id="post-filtering-for-sensible-results">3.7. Post-filtering for sensible results</h2>
<p>TSS-Predictor.pl produces a few plots to help you choose an appropriate method for cleaning up your data. The final step in the TSS-Predictor pipeline is SelectTSS.pl. This takes user-defined cutoffs and reports only transcription start sites that pass these filter(s). For SelectTSS.pl usage see the next section.</p>
<h3 id="expression-based-filtering">Expression-based filtering</h3>
<p>Transcripts that are expressed at a low level are less likely to be assembled correctly, and hence TSS predictions for the least abundant transcripts are less likely to be accurate. Expression-based filtering involves disregarding genes expressed below a user-defined level, optimal cutoffs vary between experiments and depend on a number of factors including sequencing depth and read length. In some cases expression-based filtering may not be suitable at all. TSS-Predictor uses R to produce two plots to aid the selection of appropriate cutoffs.</p>
<p>Figure 2 shows the MAD expression plot produced by TSS-Optimizer for our sample dataset, this clearly demonstrates the effect of gene expression on TSS prediction accuracy. We assume that one of the annotated transcription start sites for each expressed gene is correct. For high quality, well annotated genomes such as mouse and human this appears to be a relatively safe assumption. The MAD expression plot shows the median absolute deviation for genes, grouped into bins based on the per-base coverage (expression). We observe a much higher variation in the distances between predicted and annotated TSS for genes expressed at a lower level.</p>
<div class="figure">
<embed src="./MadExpressionPlot.pdf"><p class="caption">The MAD expression plot.</p>
</div>
<p>The dashed black line in Figure 2 shows the MAD for all genes in the dataset. TSS-Predictor advises an expression cutoff as the median per-base gene expression of the first bin below the line. The effect of trimming at this cutoff can be seen in the ExprCutoffPlot.jpg file produced by TSS-Predictor.pl, which for our example dataset is shown in Figure 3. This shows that a large number of lowly expressed outliers will be removed by disreagarding genes expressed below the red line. In principle this approach facilitates the discovery of novel TSSs by identifying highly expressed genes with predicted TSS a large distance from the nearest predicted TSS for the gene. In practise however, these 'novel' TSSs most often arise from incorrect transcript assembly, or selection of an inappropriate reference transfrag to trim. For some datasets expression-based filtering does not work well, usually because the distinction between low expression (and high variance) and high expression (and low variance) is not particuarly clear. In this case an alternative approach to filtering is distance based, see below.</p>
<div class="figure">
<embed src="./ExprCutoffPlot.pdf"><p class="caption">The expression cutoff plot. TSS-Predictor reports the value of the expression cutoff to the command line. In this instance it is 65.09, this value is printed to the terminal during TSS-Predictor analysis.</p>
</div>
<h3 id="distance-based-filtering">Distance based filtering</h3>
<p>Distance based filtering assumes explicitly that one of the annotated transcription start sites for each gene is correct. Genes with a TSS greater than a given distance from the nearest predicted TSS are disregarded. 500 bases often appears to be a good cutoff. The DistanceCutoffPlot (Figure 4) shows the effect of removing TSS greater than 500bp from the nearest annotated TSS for the gene.</p>
<div class="figure">
<embed src="./DistanceCutoffPlot.pdf"><p class="caption">The distance cutoff plot.</p>
</div>
<h2 id="selecting-a-high-confidence-tss-set-using-selecttss.pl">3.8. Selecting a high confidence TSS set using SelectTSS.pl</h2>
<p>SelectTSS.pl is used to filter TSS predictions. A tss.selected.bed file is produced which contains only transcription start sites that pass either expression-based or distance-based filtering. If both a distance and expression cutoff are supplied then only transcripts that pass both filters are reported. There is an additional option to specify whether to report the predicted TSS, or the TSS of the nearest annotated transcript for the gene.</p>
<pre><code>SelectTSS.pl [options -d -e -r] -o &lt;output prefix&gt;</code></pre>
<table>
<col width="12%" />
<col width="28%" />
<col width="58%" />
<thead>
<tr class="header">
<th align="left">Options</th>
<th align="left">Default</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">-o</td>
<td align="left">tss.selected</td>
<td align="left">Output file prefix</td>
</tr>
<tr class="even">
<td align="left">-d</td>
<td align="left">[]</td>
<td align="left">TSS distance cutoff</td>
</tr>
<tr class="odd">
<td align="left">-e</td>
<td align="left">[]</td>
<td align="left">Expression level cutoff</td>
</tr>
<tr class="even">
<td align="left">-r</td>
<td align="left">[]</td>
<td align="left">Report nearest reference transcription start site</td>
</tr>
</tbody>
</table>
<p>We want to select only high confidence TSS's in our example. The following command will report predicted TSS's that satisfy both distance AND expression criteria:-</p>
<pre><code>cd ~/...
SelectTSS.pl -d 500 -e 65.09 </code></pre>
<h1 id="output-files-and-formats-1">4.0. Output files and formats</h1>
<h2 id="coverage-file">4.1. Coverage file</h2>
<p>The coverage file produced by GetCoverage.pl is a tab-delimited text file containing the per-base coverage of each transcribed location on the reference genome as determined by Cufflinks. The file format is as follows:-</p>
<table>
<col width="20%" />
<col width="20%" />
<col width="59%" />
<thead>
<tr class="header">
<th align="left">Field position</th>
<th align="left">Field name</th>
<th align="left">Desription</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">1</td>
<td align="left">Chromosome</td>
<td align="left">Reference chromosome/scaffold</td>
</tr>
<tr class="even">
<td align="left">2</td>
<td align="left">Position</td>
<td align="left">Position on the reference chromosome</td>
</tr>
<tr class="odd">
<td align="left">3</td>
<td align="left">Read depth</td>
<td align="left">Read coverage at position on genome</td>
</tr>
</tbody>
</table>
<h2 id="trim.combined.gtf">4.2. trim.combined.gtf</h2>
<p>The trim.combined.gtf file is produced by Prep4Trim.pl. This GTF file is a reduced version of the Cuffcompare combined.gtf file and contains one transfrag record for each expressed reference gene (as determined by Cuffcompare). Specifically, the longest transfrag on the same strand as the reference gene is reported. For a detailed descripion combined.gtf format see the Cufflinks manual <a href="http://cufflinks.cbcb.umd.edu/manual.html"><code class="url">http://cufflinks.cbcb.umd.edu/manual.html</code></a>.</p>
<h2 id="trim.gtf.tmap">4.3. trim.gtf.tmap</h2>
<p>The trim.gtf.tmap file is an amended version of the transcripts.gtf.tmap file produced by Cuffcompare. This tab seperated file contains only the reference isotigs selected for trimming. There are also a couple of changes in this file relative to the Cuffcompare output. Firstly, the <em>ref_id</em> field contains the reference gene, as defined by the <em>gene_id</em> tag in the annotations GTF file, rather than the reference transcript. The second difference is the <em>gene.cov</em> field. The cuffcompare output contains a <em>cov</em> field reporting the expected transcript coverage, whereas the trim.gtf.tmap contains the expected per-base coverage for the gene calculated by summing the length-normalized coverage estimates for all transcripts in the gene (see section 5.1).</p>
<p>TSS-Predictor makes use of fields 2 and 10, the rest of the file is produced by Cufflinks. The fields in the trim.gtf tmap file are as follows:</p>
<table>
<col width="20%" />
<col width="20%" />
<col width="59%" />
<thead>
<tr class="header">
<th align="left">Field position</th>
<th align="left">Field name</th>
<th align="left">Desription</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">1</td>
<td align="left">ref_gene_id</td>
<td align="left">Gene name assigned by Cuffcompare (<em>gene_name</em> attribute of GTF file)</td>
</tr>
<tr class="even">
<td align="left">2</td>
<td align="left">ref_id</td>
<td align="left"><em>gene_id</em> attribute of the the GTF file</td>
</tr>
<tr class="odd">
<td align="left">3</td>
<td align="left">class_code</td>
<td align="left">Cuffcompare class code</td>
</tr>
<tr class="even">
<td align="left">4</td>
<td align="left">cuff_gene_id</td>
<td align="left">The gene id assigned by Cufflinks. Multiple <em>cuff_gene_id</em>'s may map to a single <em>ref_gene_id</em> attribute</td>
</tr>
<tr class="odd">
<td align="left">5</td>
<td align="left">cuff_id</td>
<td align="left">Name of the Cufflinks transcript id selected for trimming for <em>ref_gene_id</em></td>
</tr>
<tr class="even">
<td align="left">6</td>
<td align="left">FMI</td>
<td align="left">Expression of this isoform relative to major isoform for the gene (defined by Cufflnks, not used)</td>
</tr>
<tr class="odd">
<td align="left">7</td>
<td align="left">FPKM</td>
<td align="left">FPKM for the transcript <em>cuff_id</em>. Note this is the transcript (NOT gene) expression</td>
</tr>
<tr class="even">
<td align="left">8</td>
<td align="left">FPKM_conf_lo</td>
<td align="left">Upper confidence interval for FPKM</td>
</tr>
<tr class="odd">
<td align="left">9</td>
<td align="left">FPKM_conf_hi</td>
<td align="left">Lower confidence interval for FPKM</td>
</tr>
<tr class="even">
<td align="left">10</td>
<td align="left">gene.cov</td>
<td align="left">The per-base coverage estimate for the gene as used by TSS-Predictor</td>
</tr>
<tr class="odd">
<td align="left">11</td>
<td align="left">len</td>
<td align="left"><em>cuff_id</em> transcript length</td>
</tr>
<tr class="even">
<td align="left">12</td>
<td align="left">major_iso_id</td>
<td align="left"><em>cuff_id</em> of the genes (<em>cuff_gene_id</em>) major isoform</td>
</tr>
</tbody>
</table>
<h2 id="isotigs.trimmed.gtf">4.4. isotigs.trimmed.gtf</h2>
<p>The isotigs.trimmed.gtf file is an amended version of the trim.combined.gf file (see section 4.1). An additional flag, <em>contain_tss</em>, is added to field 9, this contains either &quot;yes&quot; or &quot;no&quot; and identifies the exons which are predicted to contain the primary TSS for the gene. For transcripts on the '+' strand with <em>contain_tss</em> = &quot;yes&quot; the exon start position (field 4) is amended to reflect the predicted TSS. When <em>contain_tss</em> = &quot;yes&quot; and the transcript is on the '-' strand the exon end position (field 5) is amended accordingly.</p>
<h2 id="tss.prediction.bed">4.5. tss.prediction.bed</h2>
<p>This BED (tab delimited) file contains the TSS predictions for each expressed reference gene. By default the TSS lenth is 1 base. Proximal promoter regions flanking the predicted TSS can be extracted by providing TSS-Predictor with the <em>-u</em> and <em>-d</em> flags. Please note that changing the values of <em>-u</em> and <em>-d</em> will affect the distance estimates produced in the distances.txt file. The format of tss.predictions.bed is as follows:-</p>
<table>
<col width="22%" />
<col width="22%" />
<col width="47%" />
<tbody>
<tr class="odd">
<td align="left">Field position</td>
<td align="left">Field name</td>
<td align="left">Desription</td>
</tr>
<tr class="even">
<td align="left">1</td>
<td align="left">Chromosome</td>
<td align="left">Reference chromosome/scaffold</td>
</tr>
<tr class="odd">
<td align="left">2</td>
<td align="left">Start</td>
<td align="left">TSS prediction start</td>
</tr>
<tr class="even">
<td align="left">3</td>
<td align="left">End</td>
<td align="left">TSS prediction end</td>
</tr>
<tr class="odd">
<td align="left">4</td>
<td align="left">Gene name</td>
<td align="left">The reference gene id</td>
</tr>
<tr class="even">
<td align="left">5</td>
<td align="left">Score</td>
<td align="left">Not used</td>
</tr>
<tr class="odd">
<td align="left">6</td>
<td align="left">Strand</td>
<td align="left">Strand the reference gene is on</td>
</tr>
</tbody>
</table>
<h2 id="distances.txt">4.6. distances.txt</h2>
<p>The tab-delimited distances.txt file contains a variety of information for each expressed reference gene. The format of this file is as follows:-</p>
<table>
<col width="20%" />
<col width="35%" />
<col width="44%" />
<tbody>
<tr class="odd">
<td align="left">Field position</td>
<td align="left">Field name</td>
<td align="left">Desription</td>
</tr>
<tr class="even">
<td align="left">1</td>
<td align="left">Gene</td>
<td align="left">Reference gene id</td>
</tr>
<tr class="odd">
<td align="left">2</td>
<td align="left">Chr</td>
<td align="left">Chromosome</td>
</tr>
<tr class="even">
<td align="left">3</td>
<td align="left">Prediction_start</td>
<td align="left">Predicted TSS</td>
</tr>
<tr class="odd">
<td align="left">4</td>
<td align="left">Nearest_transcript</td>
<td align="left">Closese reference transcript to TSS prediction</td>
</tr>
<tr class="even">
<td align="left">5</td>
<td align="left">Nearest_transcript_start</td>
<td align="left">TSS of nearest annotated reference transcript</td>
</tr>
<tr class="odd">
<td align="left">6</td>
<td align="left">Distance</td>
<td align="left">Distance between predicted and nearest reference TSS</td>
</tr>
</tbody>
</table>
<h2 id="tss.selected.bed">4.7. tss.selected.bed</h2>
<p>This BED file contains the subset of tss.prediction.bed containing high-confidence TSS predictions (i.e. those which have passed the filters set by SelectTSS.pl).</p>
<hr />
<h1 id="overview-of-the-algorithm-1">5.0. Overview of the algorithm</h1>
<p>TSS-Predictor works by selecting and trimming Cufflinks gene annotations based upon the ratio between the genes expected coverage (expression level) and per-base coverage at each position on the gene. Below is an overview of the steps involved:-</p>
<ol style="list-style-type: decimal">
<li><p>Firstly, reads are mapped to the reference genome allowing multi-mapping reads. Because the Cufflinks output is used by TSS-Predictor I use TopHat, but this is not essential. <strong>Depth is important for TSS prediction, for this reason I merge BAM alignments after read mapping</strong>.</p></li>
<li><p>Cufflinks is then run with the default parameters and transcribed fragments 'transfrags' are linked to reference genes using Cuffcompare. Ensembl is my favourite annotation source, but all you need is a many-to-one transcript_id to gene_id relationship in the GTF file. Note UCSC GTF files obtained from the table browser do not maintain this relationship.</p></li>
<li><p>A coverage file is produced using <em>GetCoverage.pl</em>. This script makes use of the Cufflinks output and BEDTools to calculate the per-base read coverage for every transcribed base reported in the Cufflinks cuffcompare.combined.gtf file.</p></li>
<li><p>TSS-Optimizer uses a reduced dataset of one chromosome to determine optimal conditions (see Section 3.5 <em>Using TSS-Optimizer to identify optimal coverage cutoffs</em>). TSS-Predictor predicts major TSS utilization genome-wide once the optimal conditions have been determined. TSS-Optimizer and TSS-Predictor run the same core algorithm and are described together below:-</p></li>
<li><h2>1. Prep4Trim.pl</h2></li>
</ol>
<p>Prep4Trim.pl calculates the expected coverage for each reference gene and selects one representative transfrag for trimming. Expected gene expression is calculated by summing length-normalized expression (coverage) for all Cufflinks transfrags associated with each reference gene, as defined by the <em>gene_id</em> tag in the reference GTF file (Figure 5). Additionally, one representative transfrag for each reference gene is selected for trimming. Specifically, the longest Cufflinks transfrag on the same strand as the reference gene is selected.</p>
<div class="figure">
<embed src="./Normalization.pdf"><p class="caption">Gene-wise coverage calculation. Each expressed reference gene <em>G</em> is composed of <em>n</em> transfrags (<em>T</em>), each with coverage (<em>C</em>) and length (<em>L</em>). Thus, <em>C</em><sub><em>T</em><em>i</em></sub> is the coverage for the <em>i</em><sup>th</sup> transfrag of gene <em>G</em>, and <em>L</em><sub><em>T</em><em>i</em></sub> is the length of the <em>i</em><sup>th</sup> transfrag of gene <em>G</em>.</p>
</div>
<p>Prep4Trim.pl reports two files - trim.combined.gtf and trim.gtf.tmap. Trim.combined.gtf is a reduced version of the Cufflinks cuffcompare.combined.gtf file, and contains only the transfrags selected for trimming in the next stage of the algorithm. Trim.gtf.tmap is an amended version of the cuffcompare .tmap output. This contains a mapping between each Cufflinks transfrag to be trimmed and the reference gene, as well as the amended coverage estimate for the gene (see section 4.3). All other fields in trim.gtf.tmap file are transfrag level information generated by Cuffcompare.</p>
<h2 id="trim_cufflinks_gtf.pl">5.2. Trim_Cufflinks_GTF.pl</h2>
<p>The core of the algorithm. Trim_Cufflinks_GTF.pl takes the expected gene-level coverage information from trim.gtf.tmap and transfrag information from trim.combined.gtf. Using the per-base coverage information calculated by <em>GetCoverage.pl</em> each selected isotig is trimmed at a fraction of the expected coverage for the gene (<em>X</em>).</p>
<p>The algorithm iterates, one base at a time, from the last base (3') in the first exon of each transfrag (providing it is covered by &gt;<em>X</em> * <em>C</em><sub><em>G</em></sub>) towards the first base. The GTF file is trimmed when the per-base coverage of transfrag <em>T</em><sub><em>max</em></sub> drops below <em>X</em> * <em>C</em><sub><em>G</em></sub>. If the last base of exon 1 is expressed at a lower level than <em>X</em> then the second exon is examined, so on and so forth...</p>
<p>The value of <em>X</em> is provided by the user, and determined through empirical optimization. This is what the TSS-Optimizer is for. Optimal values for <em>X</em> are then used to predict transcription start sites genome-wide using TSS-Predictor.pl.</p>
<p>Trim_Cufflinks_GTF.pl reports an isotigs.trimmed.gtf file. This is an amended version of trim.combined.gtf with revised transcription start site information. An additional variable <em>&quot;contain_tss&quot;</em> is added to the 9<sup>th</sup> field of the GTF file to identify exons predicted to contain transcription start sites.</p>
<h2 id="get-tss.pl">5.3. Get-TSS.pl</h2>
<p>This script links reports a predicted transcription start site for each reference gene in BED format. By default this is tss.prediction.bed. The <em>-u</em> and <em>-d</em> flags are used to set the number of bases upstream (<em>u</em>) and downstream (<em>d</em>) of the predicted start site to report.</p>
<h2 id="getdistances.pl">5.4. GetDistances.pl</h2>
<p>GetDistances.pl generates the distances.txt file. This tab-delimited file contains the reference gene, chromosome, predicted transcription start site, nearest reference transcript, nearest reference transcription start site and the distance between the two predictions.</p>
<h2 id="selecttss.pl">5.5. SelectTSS.pl</h2>
<p>SelectTSS.pl acts as a filter. It reporting only isoforms which are expressed above an expression cutoff defined by <em>-e</em> or a distance cutoff defined by <em>-d</em>. The distance metric used is the distance between the predicted TSS and the nearest annotated TSS for the associated gene. SelectTSS.pl makes use of the trim.gtf.tmap and distances.txt files to obtain expression and distance information respectively and reports two files tss.selected.bed amd tss.selected.distances.txt. These files are identical in format to tss.prediction.bed and distances.txt.</p>
<h2 id="collatedistances.pl-optimizer-only">5.6. CollateDistances.pl (Optimizer only)</h2>
<p>This script is run called by TSS-Optimizer and produces the CutoffDistances.txt file used to determine the optimal trimming conditions.</p>
<hr />
</body>
</html>
