---
title: "Introuducing FASTQ format and error profiles"
author: "Mary Piper, Meeta Mistry"
date: "Tuesday, November 10, 2015"
---

Approximate time: 40 minutes

## Learning Objectives:

* Setting up your project space for an NGS workflow
* Introducing bioinformatics workflows and data standards
* Learning about the FASTQ format
* Understanding quality scores
* Examining error profiles for QC analysis

## Getting your project started

Project organization is one of the most important parts of a sequencing project, but is often overlooked in the excitement to get a first look at new data. While it's best to get yourself organized before you begin analysis, it's never too late to start.

You should approach your sequencing project in a very similar way to how you do a biological experiment, and ideally, begins with experimental design. **We're going to assume that you've already designed a beautiful sequencing experiment to address your biological question, collected appropriate samples, and that you have enough statistical power.** We will cover experimental design in more detail as they apply to the different workflows. 

For all of those steps taken in the wetlab (collecting specimens, extracting DNA, prepping your samples) you've likely kept a lab notebook that details how and why you did each step, but documentation doesn't stop at the sequencer! With the computational part of your experiment there are also a ton of things you will need to keep track of. You will need to keep track of the different tools (and versions) used at every stage, in addition to the parameters you set for each. 

When working with NGS data, the raw reads you get off of the sequencer will need to pass through a number of  different tools in order to generate your final desired output. The execution of this set of tools in a specified order is commonly referred to as a *workflow* or a *pipeline*. Every step in the workflow is going to spawn many files, and inevitability, you'll need to run some of those analyses again. **Genomics projects can quickly accumulate hundreds of files across tens of folders.** Do you remember what PCR conditions you used to create your sequencing library? Probably not. Similarly, you probably won't remember whether your best alignment results were in Analysis1, AnalysisRedone, or AnalysisRedone2; or which quality cutoff you used.

Luckily, recording your computational experiments is even easier than recording wet lab data. Let's look at the best practices for organizing your genomics project. 

Your future self will thank you.


### Setting up the filesystem

To prepare for the overwhelming number of files that will be generated during your NGS workflow, it is best to set up a directory structure such that you have a designated place for files when you encounter them. In this next exercise we will setup a directory structure for the project we will be using over the next few days.

First, make sure that you are in your home directory,

```
$ pwd
```
this should give the result: `/home/user_name`

**Tip** If you were not in your home directory, the easiest way to get there is to enter the command `cd` - which always returns you to home. 

Now, make a directory for the RNA-seq analysis within the `ngs_course` folder using the `mkdir` command

```
$ mkdir ngs_course/rnaseq
```

Next you want to set up the following structure within your project directory to keep files organized:

```
rnaseq/
  ├── data/
  ├── meta/
  ├── results/
  └── logs/

```

*This is a generic structure and can be tweaked based on personal preferences.* A brief description of what might be contained within the different sub-directories is provided below:

* **`data/`**: This folder is usually reserved for any raw data files that you start with. For example, this is where the original FASTQ files (data you get from the sequencer) would reside. It is best practice to always have a copy of the data in a folder in it's raw form; as you will notice that usually a workflow is run a few times before we get it completely right.

* **`meta/`**: This folder contains any information that describes the samples you are using, which we often refer to as metadata. Usually this comes in the form of a tab-delimited/Excel file in which each row corresponds to a sample (listed using the filename for that sample in the raw data collection), and columns that follow would contain any other pertinent information for the sample (i.e sample class, demographic factors, sequencer specific information). An example of a metadata file is shown below:

![metadata](../img/metadata_example.png)

* **`results/`**: This folder will contain the output from the different tools you implement in your workflow. In some cases, you will simply have the results file as ouput but with other tools you will find a large number of intermediate files are generated. To stay organized, you should create sub-folders specific to each tool/step of the workflow. 

* **`logs/`**: It is important to keep track of the commands you run and the specific pararmeters you used, but also to have a record of any standard output that is generated while running the command. This will allow you to go back to your recorded logfiles to explore additional information (e.g., how many adapters were removed, how many reads did not align). Different tools have different ways of reporting log messages and you might have to experiment a bit to figure out what output to capture: you can redirect standard output with the `>` symbol which is equivalent to `1> (standard out)`; other tools might require you to use `2>` to re-direct the standard error instead. 
 

Let's create the directory structure by changing into `rnaseq` and then using `mkdir` to create the four directories.

```
$ cd ngs_course/rnaseq
$ mkdir data meta results logs
``` 

Verify that you have created the directories:

```
$ ls -lF
```
 
If you have created these directories, you should get the following output from that command:

```
drwxrwsr-x 2 rsk27 rsk27   0 Jun 17 11:21 data/
drwxrwsr-x 2 rsk27 rsk27   0 Jun 17 11:21 logs/
drwxrwsr-x 2 rsk27 rsk27   0 Jun 17 11:21 meta/
drwxrwsr-x 2 rsk27 rsk27   0 Jun 17 11:21 results/

```
Now we will create the subdirectories to setup for our RNA-Seq analysis, and populate them with data where we can. The first step will be checking the quality of our data, and trimming the files if necessary. We need to create two directories within the `data` directory, one folder for untrimmed reads and another for our trimmed reads: 

```
$ cd ~/ngs_course/rnaseq/data
$ mkdir untrimmed_fastq
$ mkdir trimmed_fastq
```
    
The raw_fastq data we will be working with is currently in the `unix_lesson/raw_fastq` directory. We need to copy the raw fastq files to our `untrimmed_fastq` directory:

`$ cp ~/ngs_course/unix_lesson/raw_fastq/*fq  ~/ngs_course/rnaseq/data/untrimmed_fastq`

Later in the workflow when we perform alignment, we will require reference files to map against. These files are also in the `unix_lesson` directory, you can copy the entire folder over into `data`:

`$ cp -r ~/ngs_course/unix_lesson/reference_data rnaseq/data`
 
OK we're all setup for our analysis!

## Bioinformatics workflows and data standards

An example of the workflow we will be using for our RNA-seq analysis is provided below with a brief description of each step. As we work through this pipeline during Sessions I and II, we will discuss each in more detail.


1. Quality control - Assessing quality using FastQC
2. Quality control - Trimming and/or filtering reads (if necessary)
3. Align reads to reference genome using STAR (splice-aware aligner)
4. Count the number of reads mapping to each gene using htseq-count
5. Statistical analysis (count normalization, linear modeling using R-based tools)

![workflow](../img/rnaseq_workflow.png)

These workflows in bioinformatics adopt a plug-and-play approach in that the output of one tool can be easily used as input to another tool without any extensive configuration. Having standards for data formats is what makes this feasible. Standards ensure that data is stored in a way that is generally accepted and agreed upon within the community. The tools that are used to analyze data at different stages of the workflow are therefore built under the assumption that the data will be provided in a specific format. We will discuss the various data formats as we encounter them within the RNA-seq workflow.


### FASTA format

We have all encountered a FASTA file format at some point during our research career. FASTA was first introduced in the late 80's and has now become a standard in bioinformatics. It is a **text-based file format used to represent either protein or nucleotide sequences.** Each record in a FASTA file is represented by a `>` character, as this marks the header. The header can contain any information about the sequence (i.e. accession, sequence name, description). Any lines that follow the header contain sequence information.

```
>gi|340780744|ref|NC_015850.1| Acidithiobacillus caldus SM-1 chromosome, complete genome 
ATGAGTAGTCATTCAGCGCCGACAGCGTTGCAAGATGGAGCCGCGCTGTGGTCCGCCCTATGCGTCCAACTGGAGCTCGTCACGAGTCCGCAGCAGTT
CAATACCTGGCTGCGGCCCCTGCGTGGCGAATTGCAGGGTCATGAGCTGCGCCTGCTCGCCCCCAATCCCTTCGTCCGCGACTGGGTGCGTGAACGCA
TGGCCGAACTCGTCAAGGAACAGCTGCAGCGGATCGCTCCGGGTTTTGAGCTGGTCTTCGCTCTGGACGAAGAGGCAGCAGCGGCGACATCGGCACCG
ACCGCGAGCATTGCGCCCGAGCGCAGCAGCGCACCCGGTGGTCACCGCCTCAACCCAGCCTTCAACTTCCAGTCCTACGTCGAAGGGAAGTCCAATCA
GCTCGCCCTGGCGGCAGCCCGCCAGGTTGCCCAGCATCCAGGCAAATCCTACAACCCACTGTACATTTATGGTGGTGTGGGCCTCGGCAAGACGCACC

>gi|129295|sp|P01013|OVAX_CHICK GENE X PROTEIN (OVALBUMIN-RELATED)
QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE

```

### FASTQ format: FASTA + Quality scores

The [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file format is the defacto file format for sequence reads generated from next-generation sequencing technologies. This file format evolved from FASTA in that it contains sequence data, but also contains quality information. Similar to FASTA, the FASTQ file begins with a header line. The difference is that the FASTQ header is denoted by a `@` character. For a single record (sequence read) there are four lines, each of which are described below:

|Line|Description|
|----|-----------|
|1|Always begins with '@' and then information about the read|
|2|The actual DNA sequence|
|3|Always begins with a '+' and sometimes the same info in line 1|
|4|Has a string of characters which represent the quality scores; must have same number of characters as line 2|

Let's use an example read from our data set:

```
@HWI-ST330:304:H045HADXX:1:1101:1111:61397
CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
+
@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################
```

As mentioned previously, line 4 has characters encoding the quality of each nucleotide in the read. The legend below provides the mapping of quality scores (Phred-33) to the quality encoding characters. ** *Different quality encoding scales exist (differing by offset in the ASCII table), but note the most commonly used one is fastqsanger* **

 ```
 Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI
                   |         |         |         |         |
    Quality score: 0........10........20........30........40                                
```
 
Using the quality encoding character legend, the first nucelotide in the read (C) is called with a quality score of 31 and our Ns are called with a score of 2. **As you can tell by now, this is one of our bad reads.** 

Each quality score represents the probability that the corresponding nucleotide call is incorrect. This quality score is logarithmically based and is calculated as:

	Q = -10 x log10(P), where P is the probability that a base call is erroneous

These probabaility values are the results from the base calling algorithm and dependent on how much signal was captured for the base incorporation. The score values can be interpreted as follows:

|Phred Quality Score |Probability of incorrect base call |Base call accuracy|
|:-------------------|:---------------------------------:|-----------------:|
|10	|1 in 10 |	90%|
|20	|1 in 100|	99%|
|30	|1 in 1000|	99.9%|
|40	|1 in 10,000|	99.99%|
|50	|1 in 100,000|	99.999%|
|60	|1 in 1,000,000|	99.9999%|

Therefore, for the first nucleotide in the read (C), there is less than a 1 in 1000 chance that the base was called incorrectly. Whereas, for the the end of the read there is greater than 50% probabaility that the base is called incorrectly.

#### Assessing quality

The quality scores are useful in determining whether a sample is good or bad. Rather than looking at quality scores for each individual read, we use a tool called [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to looks at quality collectively across all reads within a sample. The image below is a plot that indicates a (very) good quality sample:

![good_quality](../img/good_quality.png)

On the x-axis you have the base position in the read, and on the y-axis you have quality scores. In this example, the sample contains reads that are 40 bp long. For each position, there is a box plotted to illustrate the distribution of values (with the whiskers indicating the lowest value observed at that position). For every position here, the quality values do not drop much lower than 32 -- which if you refer to the table above is a pretty good quality score. The plot background is also color-coded to identify good (green), acceptable (yellow), and bad (red) quality scores.  


Now let's take a look at a quality plot on the other end of the spectrum. 

![good_quality](../img/bad_quality.png)

Here, we see positions within the read in which the boxes span a much wider range. Also, quality scores drop quite low into the 'bad' range, particularly on the tail end of the reads. When you encounter a quality plot such as this one, the first step is to troubleshoot. Why might we be seeing something like this? 

#### Error profiles

The FASTQC tool produces several other diagostic plots to assess sample quality, in addition to the one plotted above. We will discuss various error profiles and potential causes in the [slides provided here](../slides/error_profiles_mm.pdf). 

----

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
