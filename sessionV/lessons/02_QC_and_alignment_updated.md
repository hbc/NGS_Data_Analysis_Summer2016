---
title: "ChIP-Seq QC and alignment"
author: "Mary Piper, Radhika Khetani"
date: "Monday, April 4th, 2016"
---

Contributors: Mary Piper, Radhika Khetani

Approximate time: 1 hour

## Learning Objectives

* to use previous knowledge of quality control steps to perform FastQC and trimming
* to understand parameters and perform alignment using Bowtie2

# ChIP-Seq analysis 

Now that we have our files and directory structure, we are ready to begin our ChIP-Seq analysis. 

![workflow_QC](../img/chipseq_workflow_QC_partial.png)

## Quality Control
For any NGS analysis method, our first step is ensuring our reads are of good quality prior to aligning them to the reference genome. We will use FastQC to get a good idea of the overall quality of our data, to identify whether any samples appear to be outliers, to examine our data for contamination, and to determine a trimming strategy.

Let's run FastQC on all of our files. 

Start an interactive session if you are not already in one and change directories to the untrimmed_fastq folder.

```
$ cd ~/ngs_course/chipseq/data/untrimmed_fastq 

$ module load seq/fastqc/0.11.3 

$ fastqc H1hesc_Input_Rep1_chr12.fastq 
```

Now, move all of the `fastqc` files to the `results/untrimmed_fastqc` directory:

`$ mv *fastqc* ../../results/untrimmed_fastqc/`

Transfer the FastQC zip file for Input replicate 1 to your local machine using FileZilla and view the report.

![fastqc](../img/fastqc_input_rep1.png)

Based on the sequence quality plot, trimming should be performed from both ends of the sequences. We will use Trimmomatic to trim the reads from both ends of the sequence using the following parameters:

* `SE`: Single End reads
* `-threads`: number of threads / cores
* `-phread33`: quality score format
* `LEADING`: cut bases off the start of a read, if below a threshold quality
* `TRAILING`: cut bases off the end of a read, if below a threshold quality
* `MINLEN`: drop an entire read if it is below a specified length

Since we are only trimming a single file, we will run the command in the interactive session rather than creating a script:

```
$ java -jar /opt/Trimmomatic-0.33/trimmomatic-0.33.jar SE \
-threads 4 \
-phred33 \
H1hesc_Input_Rep1_chr12.fastq \
../trimmed_fastq/H1hesc_Input_Rep1_chr12.qualtrim20.minlen36.fq \
LEADING:20 \
TRAILING:20 \
MINLEN:36
```

Let's see how much trimming improved our reads by running FastQC again:

`$ fastqc ../trimmed_fastq/H1hesc_Input_Rep1_chr12.qualtrim20.minlen36.fq`

Move the FastQC folders to the results directory for trimmed FastQC results:

`$ mv ../trimmed_fastq/*fastqc* ../../results/trimmed_fastqc/`

Using Filezilla, transfer the file for the trimmed Input replicate 1 FastQC to your computer.

![trimmed_fastqc](../img/chipseq_trimmed_fastqc.png)

## Alignment

![workflow_align](../img/chipseq_workflow_align_partial.png)

Now that we have removed the poor quality sequences from our data, we are ready to align the reads to the reference genome. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) is a fast and accurate alignment tool that indexes the genome with an FM Index based on the Burrows-Wheeler Transform to keep memory requirements low for the alignment process. 

*Bowtie2* supports gapped, local and paired-end alignment modes. It works best for reads that are at least 50 bp (shorter read lengths should use Bowtie1), and it can perform soft-clipping to remove poor quality bases [[1](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2009-10-3-r25)].

To perform peak calling for ChIP-Seq analysis, we need our alignment files to contain only **uniquely mapping reads** (no multi-mappers or duplicate reads) in order to increase confidence in site discovery and improve reproducibility. Since there is no parameter in Bowtie2 to keep only uniquely mapping reads, we will need to perform the following steps to generate alignment files containing only the uniquely mapping reads:

1. Create a bowtie2 index
2. Align reads with bowtie2 and output a SAM alignment file
3. Change alignment file format from SAM to BAM
4. Sort BAM file by read coordinate locations
5. Filter to keep only uniquely mapping reads (this will also remove any unmapped reads)

### Creating Bowtie2 index

To perform the Bowtie2 alignment, a genome index is required. **We previously generated the genome indexes for you**, and they exist in the `reference_data` directory.

However, if you needed to create a genome index yourself, you would use the following command:

```
# DO NOT RUN

bowtie2-build <path_to_reference_genome.fa> <prefix_to_name_indexes>

```

### Aligning reads with Bowtie2

Since we have our indexes already created, we can get started with read alignment. Change directories to the `bowtie2` folder:

```

$ cd ~/ngs_course/chipseq/results/bowtie2

```

We will perform alignment on our single trimmed sample, `H1hesc_Input_Rep1_chr12.qualtrim20.minlen36.fq`. Details on Bowtie2 and its functionality can be found in the [user manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml); we encourage you to peruse through to get familiar with all available options.

The basic options for aligning reads to the genome using Bowtie2 are:

* `-p`: number of processors / cores
* `-q`: reads are in FASTQ format
* `-x`: /path/to/genome_indices_directory
* `-U`: /path/to/FASTQ_file
* `-S`: /path/to/output/SAM_file

```
$ bowtie2 -p 6 -q \
-x ~/ngs_course/chipseq/data/reference_data/chr12 \
-U ~/ngs_course/chipseq/data/trimmed_fastq/H1hesc_Input_Rep1_chr12.qualtrim20.minlen36.fq \
-S ~/ngs_course/chipseq/results/bowtie2/H1hesc_Input_Rep1_chr12_aln_unsorted.sam

```
>*NOTE: If you had added the bcbio path ino your `.bashrc` file last session you should be able to use bowtie2 it without loading a module. If not, load the module using `module load seq/bowtie/2.2.4`
>

### Changing file format from SAM to BAM

While the SAM alignment file output by Bowtie2 is human readable, we need a BAM alignment file for downstream tools. Therefore, we will use [Samtools](http://samtools.github.io) to convert the file formats. The command we will use is `samtools view` with the following parameters

* `-h`: include header in output
* `-S`: input is in SAM format
* `-b`: output BAM format
* `-o`: /path/to/output/file

```
$ samtools view -h -S \
-b H1hesc_Input_Rep1_chr12_aln_unsorted.sam \
-o H1hesc_Input_Rep1_chr12_aln_unsorted.bam
```

### Sorting BAM files by genomic coordinates

Before we can filter to keep the uniquely mapping reads, we need to sort our BAM alignment files by genomic coordinates. To perform this sort, we will use [Sambamba](http://lomereiter.github.io/sambamba/index.html), which is a tool that quickly processes BAM and SAM files. It is similar to SAMtools, but has unique functionality.
The command we will use is `sambamba sort` with the following parameters:

* `-t`: number of threads / cores
* `-o`: /path/to/output/file

```
$ sambamba sort -t 6 \
-o H1hesc_Input_Rep1_chr12_aln_sorted.bam \
H1hesc_Input_Rep1_chr12_aln_unsorted.bam 
```

### Filtering uniquely mapping reads

Finally, we can filter the uniquely mapped reads. We will use the `sambamba view` command with the following parameters:

* `-t`: number of threads / cores
* `-h`: print SAM header before reads
* `-f`: format of output file (default is SAM)
* `-F`: set [custom filter](https://github.com/lomereiter/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax)

```
$ sambamba view -h -t 6 -f bam \
-F "[XS] == null and not unmapped " H1hesc_Input_Rep1_chr12_aln_sorted.bam > H1hesc_Input_Rep1_chr12_aln.bam
```

Now that the alignment files contain only uniquely mapping reads, we can assess the quality of our alignment for peak calling.

## Alignment quality assessment

### *phantompeakqualtools* for analysis of the cross-correlation peak and relative phantom peak

The *[phantompeakqualtools](https://code.google.com/archive/p/phantompeakqualtools/)* package allows for the generation of enrichment and quality measures for ChIP-Seq data [[1](http://www.g3journal.org/content/4/2/209.full)]. We will be using the package to compute the predominant insert-size (fragment length) based on strand cross-correlation peak and data quality measures based on relative phantom peak.

#### Set up

The package is written as an R script, `run_spp.R`, with various options that can be specified when running from the command line. It is run using `R` and `samtools`. To get set up we will need to start an interactive session, load the necessary modules and set up the directory structure:

```
$ bsub -Is -n 4 -q interactive bash

$ module load stats/R/3.2.1 seq/samtools/1.2

$ cd ~/ngs_course/chipseq/results/bowtie2

$ mkdir phantompeaks

$ cd phantompeaks
```

We will need to install the R library, `caTools` to run the script, so let's [set up a personal R library on Orchestra](https://wiki.med.harvard.edu/Orchestra/PersonalRPackages):

```
$ mkdir -p ~/R/library

$ echo 'R_LIBS_USER="~/R/library"' >  $HOME/.Renviron

$ export R_LIBS_USER="/home/user123/R/library"
```

Now that we have our personal R library set up, we can install packages using R:

```
$ R
```

In R, use the install.packages() function to install `caTools`:

```
> install.packages("caTools", lib="~/R/library")

> quit()

```

#### Downloading *phantompeakqualtools*

To use this *phantompeakqualtools* package, we need to download it from the project website. On the [project website](https://code.google.com/archive/p/phantompeakqualtools/), click on the *Downloads* option on the left-hand side of the page. The *Downloads* page has all updates for the package, with the most recent being from 2013. 

Right-click on the link for the most recent update, and copy the link.

Download the *phantompeakqualtools* to your directory:

```
$ wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/phantompeakqualtools/ccQualityControl.v.1.1.tar.gz

$ ls
```

You should see `ccQualityControl.v.1.1.tar.gz` appear in the folder. This is a compressed folder, to extract the contents we use the `tar -xzf` command:

```
$ tar -xzf ccQualityControl.v.1.1.tar.gz
```

You should now see a `phantompeakqualtools` folder. Let's explore the contents a bit:

```
$ cd phantompeakqualtools

$ ls -l
```

Note the script for generating the quality metrics, `run_spp.R`. Let's check out the `README.txt`:

```
$ less README.txt
```

Within the README.txt are all the commands, options, and output descriptions.

#### Running *phantompeakqualtools*

To determine strand cross-correlation peak / predominant fragment length OR print out quality measures, typical usage for running the `run_spp.R` script from the command line include:

* `-c`: full path and name (or URL) of tagAlign/BAM file
* `-savp`: save cross-correlation plot
* `-out`: will create and/or append to a file named several important characteristics of the dataset described in more detail below.

```
$ Rscript run_spp.R -c=<tagAlign/BAMfile> -savp -out=<outFile>
```

From the `phantompeakqualtools` directory, run a 'for loop' to run the script on every BAM file:

```
$ for bam in ../../*aln.bam
do 
Rscript run_spp.R -c=$bam -savp -out=${bam}.qual >${bam}.Rout
done
```

#### Output for *phantompeakqualtools*

Let's explore the output files:

##### .Rout and .qual files

The *.Rout and *.qual files contain quality measure information
```
$ less ../../H1hesc_Nanog_Rep1_chr12_aln.bam.Rout
$ less ../../H1hesc_Nanog_Rep1_chr12_aln.bam.qual
```

The output file is tab-delimited with the columns containing the following information:

- COL1: Filename: tagAlign/BAM filename 
- COL2: numReads: effective sequencing depth i.e. total number of mapped reads in input file 
- COL3: estFragLen: comma separated strand cross-correlation peak(s) in decreasing order of correlation. (**NOTE:** The top 3 local maxima locations that are within 90% of the maximum cross-correlation value are output. In almost all cases, the top (first) value in the list represents the predominant fragment length.) 
- COL4: corr_estFragLen: comma separated strand cross-correlation value(s) in decreasing order (col2 follows the same order) 
- COL5: phantomPeak: Read length/phantom peak strand shift 
- COL6: corr_phantomPeak: Correlation value at phantom peak 
- COL7: argmin_corr: strand shift at which cross-correlation is lowest 
- COL8: min_corr: minimum value of cross-correlation 
- COL9: Normalized strand cross-correlation coefficient (NSC) = COL4 / COL8 
- COL10: Relative strand cross-correlation coefficient (RSC) = (COL4 - COL8) / (COL6 - COL8) 
- COL11: QualityTag: Quality tag based on thresholded RSC (codes: -2:veryLow,-1:Low,0:Medium,1:High,2:veryHigh)

Two of the more important values to observe are the NSC and RSC values:

**NSC:** values range from a minimum of 1 to larger positive numbers. 1.1 is the critical threshold. Datasets with NSC values much less than 1.1 (< 1.05) tend to have low signal to noise or few peaks (this could be biological eg.a factor that truly binds only a few sites in a particular tissue type OR it could be due to poor quality)

**RSC:** values range from 0 to larger positive values. 1 is the critical threshold. RSC values significantly lower than 1 (< 0.8) tend to have low signal to noise. The low scores can be due to failed and poor quality ChIP, low read sequence quality and hence lots of mismappings, shallow sequencing depth (significantly below saturation) or a combination of these. Like the NSC, datasets with few binding sites (< 200) which is biologically justifiable also show low RSC scores.

##### .pdf files
The cross-correlation plots show the best estimate for strand shift and the cross-correlation values. This file can be viewed by transferring it to your local machine using FileZilla. Copy `H1hesc_Nanog_Rep1_chr12_aln.pdf` to your machine to view the strand shift.

### *deepTools* for quality assessment of read coverage

Using the *[deepTools](http://deeptools.readthedocs.org/en/latest/content/list_of_tools.html)*, suite of tools, we can assess the quality of our alignments for each of our samples using several metrics.

Assessing and visualizing alignment quality using *deepTools* requires three steps: 

1. Indexing the BAM alignment files
2. Calculation of the read coverage scores using the `multiBamSummary` tool
3. Visualizing how read coverage scores compare between samples

#### Indexing the BAM alignment files

Similar to this step in previous lessons, we will index our BAM files using the `samtools index` tool.

Since we loaded `samtools` to use the *phantompeakqualtools*, we do not need to load it again. Let's just create a `deeptools` directory within the `bowtie2` folder:

```
$ cd ~/ngs_course/chipseq/results/bowtie2

$ mkdir deeptools 

$ cd deeptools

```

Then, we can index the BAM files by using the command: `samtools index path/to/bam`. 

Since we would like to index all of our BAM files containing uniquely mapping reads in the `bowtie2` folder, we can use a 'for loop' to index all files ending with `aln.bam`:

```
$ for bam in ../*aln.bam
> do
> samtools index $bam
> done
```
Now we should have an index (BAI) file for each of our BAM files.

#### Calculation of the read coverage scores using the `multiBamSummary` tool

The `multiBamSummary` tool will calculate the read coverage scores for specific genomic regions between samples and provide the output as a binary compressed numpy array (.npz) file. Alternatively, the analysis can be performed on the entire genome by changing the mode of this tool to ‘bins’.

```
multiBamSummary bins --ignoreDuplicates -p 6 \
--bamfiles *aln.bam \
-out deeptools_multiBAM.out.npz \
--outRawCounts readCounts.tab
```

#### Visualizing how read coverage quality metrics

Now that we have the read coverage scores calculated for all samples, we can now analyze the coverage between samples using a variety of the *deepTools* tools:

##### 1. Sample correlation - `plotCorrelation` tool

The `plotCorrelation` tool allows us to visualize the similarity between samples based on their read coverage of regions of the genome. For example, we can compare two samples to determine whether they have similar coverage profiles with either a heatmap or a scatterplot:

![correlate](../img/QC_bamCorrelate_deeptools.png)

We can analyze read coverage similarity using heatmap to perform heirarchical clustering and determine whether our sample groups cluster well (i.e. have similar read coverage profiles within and between sample groups).

```
plotCorrelation --corData deeptools_multiBAM.out.npz \
--plotFile deeptools_heatmap.png \
--corMethod pearson \
--whatToPlot heatmap \
--labels [Input_Rep1 Input_Rep2 Nanog_Rep1 Nanog_Rep2 Pou5f1_Rep1 Pou5f1_Rep2] \
--plotNumbers
```

The same `plotCorrelation` tool can be used to examine the coverage profiles by scatterplot:

plotCorrelation --corData deeptools_multiBAM.out.npz \
--plotFile deepTools_scatterplot.png \
--corMethod pearson \
--whatToPlot scatterplot \
--labels [Input_Rep1 Input_Rep2 Nanog_Rep1 Nanog_Rep2 Pou5f1_Rep1 Pou5f1_Rep2]

##### 2. Sample variability - `plotPCA` tool

The next quality metric we will explore is the principal component analysis (PCA) of our read coverage calculations. PCA can be used to determine whether samples display greater variability between experimental conditions than between replicates of the same treatment based on information (read coverage values) from thousands of regions. PCA is also useful to identify unexpected patterns, such as those caused by batch effects or outliers. 

You will use the tool `plotPCA` to sort the principal components according to the amount of variability of the data that they explain and generate two plots:

- the PCA plot for the top two principal components eigenvalues 
- the Scree plot for the top five principal components where the bars represent the amount of variability explained by the individual factors and the red line traces the amount of variability is explained by the individual components in a cumulative manner [[1]](http://deeptools.readthedocs.org/en/latest/content/tools/plotPCA.html)

![PCA](../img/PCA_deeptools.png)

```
plotPCA --corData deeptools_multiBAM.out.npz \
--plotFile deepTools_pcaplot.png \
-T "PCA of read counts" \
--outFileNameData pcaProfile.tab \
--labels [Input_Rep1 Input_Rep2 Nanog_Rep1 Nanog_Rep2 Pou5f1_Rep1 Pou5f1_Rep2]
```
##### 3. Sample sequencing depth - `plotCoverage` tool

The `plotCoverage` tool will generate plots to explore the average number of reads per base pair in the genome. The tool will generate two plots, giving the frequencies of read coverage and the fraction of bases versus read coverage.

![coverage](../img/plotCoverage_deeptools.png)

```
plotCoverage --bamfiles *aln.bam \
--ignoreDuplicates \
-o deepTools_coverageplots.png \
--labels Input_Rep1 Input_Rep2 Nanog_Rep1 Nanog_Rep2 Pou5f1_Rep1 Pou5f1_Rep2
```
##### 4. Sample signal strength - `plotFingerprints` tool

The `plotFingerprints` tool "determines how well the signal in the ChIP-seq sample can be differentiated from the background distribution of reads in the control sample" [[2](http://deeptools.readthedocs.org/en/latest/content/tools/plotFingerprint.html)].  

"For factors that will enrich well-defined, rather narrow regions (e.g. transcription factors such as p300), the resulting plot can be used to assess the strength of a ChIP, but the broader the enrichments are to be expected, the less clear the plot will be" [[2](http://deeptools.readthedocs.org/en/latest/content/tools/plotFingerprint.html)].

The tool will generate a plot for the cumulative read coverages for each sample.

![fingerprint](../img/plotFingerprint_deeptools.png)


***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

