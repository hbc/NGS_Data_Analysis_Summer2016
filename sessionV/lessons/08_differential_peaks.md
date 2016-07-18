---
title: "Differenial Peak calling MACS2"
author: "Meeta Mistry"
date: "Sunday, March 6th, 2016"
---

Contributors: Meeta Mistry, 

Approximate time: 90 minutes

## Learning Objectives

* Learning how to use MACS2 to compare peaks between two sample groups


## Differential peak enrichment 

Increasing number of ChIP-seq experiments are investigating transcription factor binding under multiple experimental conditions, for example, various treatment conditions, several distinct time points and different treatment dosage levels. Hence, identifying differential binding sites across multiple conditions is of practical importance in biological and medical research. 

There are various methods/tools available when investigating narrow peaks, and the choice of tool will depend heavily on your experimental design. 

![diffbind](../img/diff-peaks.png)

In our case, we are interested in identifying differences in binding between two transcription factors. For each group we have two replicates, and it would be best to use tools that make use of these replicates (i.e [DiffBind](http://bioconductor.org/packages/release/bioc/html/DiffBind.html), [ChIPComp](https://www.bioconductor.org/packages/3.3/bioc/html/ChIPComp.html)) to compute statistics reflecting how significant the changes are. 


## DiffBind

DiffBind is an R package that is used for identifying sites that are differentially bound between two sample groups. It works primarily with sets of peak calls ('peaksets'), which are sets of genomic intervals representing candidate protein binding sites for each sample. It includes functions that support the processing of peaksets, including overlapping and merging peak sets across an entire dataset, counting sequencing reads in overlapping intervals in peak sets, and identifying statistically significantly differentially bound sites based on evidence of binding affinity (measured by differences in read densities). We will discuss the importance of each step but for more information take a look at the [DiffBind vignette.](http://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf).


### Setting up

#### X11 forwarding

A number of programs with graphical user interfaces use the X11 system which lets the program **run on an Orchestra computer, but show the graphics on your desktop**. R is one of these programs. For this lesson we are going to want to plot diagnostic figures as we walk through the workflow, and be able to look at them interactively (instead of saving each to file and view locally).

To do this, you need to have an **X11 server running on your desktop**, and your SSH connection needs to **have X11 forwarding enabled on Orchestra**.

> *NOTE:* the X11 server is the program on your laptop that drives the user's display and handles connections from X11 clients. If you are using a Mac, this will be [XQuartz](https://www.xquartz.org/) and for PC users this would be [Xming](https://sourceforge.net/projects/xming/). You should already have these installed on your laptops. 

To setup X11 forwarding on Orchestra we need to list the settings in your SSH client's configuration file. Login to Orchestra and open up the config file using `vim`:

	$ ssh ecommons_id@orchestra.med.harvard.edu
	$ vim ~/.ssh/config
	
Now type in the following and save an exit:

	Host orchestra.med.harvard.edu
	ForwardX11 Yes


You're all setup! Log back in to Orchestra using the `-X` parameter to enable X11 forwarding:


	$ ssh -X ecommons_id@orchestra.med.harvard.edu


Then start up an interactive session with 4 cores:

	$ bsub -Is -n 4 -q interactive bash

Let's load the R module. We are going to use version 3.2.1 since it has DiffBind installed for us. 

	$ module load stats/R/3.2.1

Navigate to the `results` directory we have been working in and create a new directory for our DiffBind analysis:

	$ cd ~/ngs_course/chipseq/results
	$ mkdir diffBind


Finally, you will need the **sample sheet** which contains metadata information. Copy this over to your `diffBind` directory and then we will take a quick look at what is contained in it.

	$ cp /groups/hbctraining/ngs-data-analysisSummer2016/chipseq/ENCODE/diffBind/samples_DiffBind.csv diffBind/

	$ less diffBind/samples_DiffBind.csv


The **sample sheet** contains a row for each peak set (which in most cases is every ChIP sample) and several columns of required information, which allows us to easily load the associated data in one single command. _The column headers have specific names that are expected by DiffBind_. 

* SampleID: Identifier string for sample
* Replicate: Replicate number of sample
* Tissue, Factor, Condition, Treatment: Identifier strings for up to four different factors (need one of these at minimum)
* bamReads: file path for bam file containing aligned reads for ChIP sample
* bamControl: file path for bam file containing aligned reads for control sample
* ControlID: Identifier string for control sample
* Peaks: path for file containing peaks for sample
* PeakCaller: Identifier string for peak caller used. Possible values include “raw”, “bed”, “narrow”, “macs”


> _NOTE:_ The paths provided in the sample sheet for the alignment files (BAM) and peak calls (narrowPeak) were generated using the full dataset and are pointing to our shared directory. In this way we do not have to copy over the files. BAM files were downloaded directly from ENCODE. Links to the full dataset files are provided in the [QC markdown](https://github.com/hbc/NGS_Data_Analysis_Summer2016/blob/master/sessionV/lessons/07_IDR_assessing_replicates.md#running-idr). Peaks were called using MACS2 default parameters with `-q  0.05`.


Finally, let's open up R and load the required libraries:

	$ R

```
R version 3.2.1 (2015-06-18) -- "World-Famous Astronaut"
Copyright (C) 2015 The R Foundation for Statistical Computing
Platform: x86_64-unknown-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> 

```

	> library(DiffBind)
	



### Reading in Peaksets

The first step is to read in a set of peaksets and associated metadata. This is done using the sample sheet. Once the peaksets are read in, a merging function finds all overlapping peaks and derives a single set of unique genomic intervals covering all the supplied peaks (a consensus peakset for the experiment). *A region is considered for the consensus set if it appears in more than two of the samples.*

```
samples <- read.csv('diffBind/samples_DiffBind.csv')
dbObj <- dba(sampleSheet=samples)

```

Take a look at what information gets summarized in the `dbObj`. *How many consensus sites were identified for this dataset? Which sample has a disproportionatley larger number of peaks?*

	dbObj
	
	4 Samples, 5558 sites in matrix (17547 total):
           ID Factor Replicate Caller Intervals
	1  Nanog-Rep1  Nanog         1 narrow     11047
	2  Nanog-Rep2  Nanog         2 narrow      3562
	3 Pou5f1-Rep1 Pou5f1         1 narrow      8995	4 Pou5f1-Rep2 Pou5f1         2 narrow      3562

### Occupancy analysis:


***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*







