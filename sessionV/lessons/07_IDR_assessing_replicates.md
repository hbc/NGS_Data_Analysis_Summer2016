---
title: "Assessing reproducibility between replicates"
author: "Meeta Mistry"
date: "Monday, June 13th, 2016"
---

Contributors: Meeta Mistry, 

Approximate time: 90 minutes

## Learning Objectives

* Learning how to use IDR


## Handling replicates in ChIP-Seq
 
As with any high-throughput experiment, any single assay is often subject to a substantial amount of variability. Thus, it is highly recommended to setup your experimental design with a minimum of 2-3 biological replicates. Presumably, two replicates measuring the same underlying biology should have high consistency but that is not always the case. In order to evaluate consistency between replicates we require metrics that objectively assess the reproducibility of high-throughput assays.

In our case, we have two replicates for each transcription factor. We want to consider the peaks that are consistent in both replicates before we can compare the peaks from the two transcription factors to one another.

<img src=../img/idr_samples.png width=800> 


Common methods for handling replicates includes taking overlapping peak calls across replicates and then assessing differences in binding regions. However, these are simple methods that do not employ any statistical testing and so we know little about how robust these peaks truly are.

> **Historical Note:** A simpler heuristic for establishing reproducibility was previously used as a standard for depositing ENCODE data and was in effect when much of the currently available data was submitted. According to this standard, either 80% of the top 40% of the targets identified from one replicate using an acceptable scoring method should overlap the list of targets from the other replicate, or target lists scored using all available reads from each replicate should share more than 75% of targets in common. As with the current standards, this was developed based on experience with accumulated ENCODE ChIP-seq data, albeit with a much smaller sample size.


## Irreproducibility Discovery Rate (IDR)

[IDR](https://sites.google.com/site/anshulkundaje/projects/idr) is a framework developed by Qunhua Li and Peter Bickel's group that **compares a pair of ranked lists of regions/peaks and assigns values that reflect its reproducibility.** 

<img src=../img/idr_figure.png> 

It is extensively used by the ENCODE and modENCODE projects and is part of their [ChIP-seq guidelines and standards](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/). It has been established for submission of ChIP-seq data sets and have been constructed based on the historical experiences of ENCODE ChIP-seq data production groups.

### Why IDR?

* IDR avoids choices of initial cutoffs, which are not comparable for different callers 
* IDR does not depend on arbitrary thresholds and so all regions/peaks are considered. 
* It is based on ranks, so does not require the input signals to be calibrated or with a specific fixed scale (only order matters).


### Components of IDR

The IDR approach creates a curve, from which it then quantitatively assesses when the ﬁndings are no longer consistent across replicates. There are three main components: 

1) A **correspondence curve**: a graphical representation of matched peaks as you go down the ranked list. Qualitative, not adequate for selecting signals.

<img src=../img/corr_curve.png width=400> 

2) An **inference procedure**: summarizes the proportion of reproducible and irreproducible signals. Quantitative, using a copula mixture model.

> What proportion of identifications have a poor correspondence, i.e. falling into ”noise”?
> How consistent are the identifications before reaching breakdown?

3) **Irreproducible Discovery Rate (IDR)**: Derive a significance value from the inference procedure (#2) in a fashion similar to FDR, and can be used to control the level of irreproducibility rate when selecting signals.
i.e. 0.05 IDR means that peak has a 5% chance of being an irreproducible discovery



### The IDR pipeline

There are three main steps to the IDR pipeline:

1. Evaluate peak consistency between true replicates
2. Evaluate peak consistency between pooled pseudo-replicates
3. Evaluate self-consistency for each individual replicate

<img src=../img/idr_pipeline.png> 

> This figures is taken from the [ENCODE ChIP-Seq Guidelines](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/#box3).

_We will only be running Step 1 in this lesson, but will discuss steps 2 and 3 in a bit more detail._


## Running IDR

To run IDR we will need to use the full dataset. We started with BAM files downloaded from ENCODE using the links provided below:

* Input: https://www.encodeproject.org/experiments/ENCSR000BHL/
* Nanog Replicates: https://www.encodeproject.org/experiments/ENCSR000BMT/
* Pou5f1 Replicates: https://www.encodeproject.org/experiments/ENCSR000BMU/


Using MACS2, we **called peaks using slightly looser thresholds (p < 0.001)** than we would normally use for peak calling. This is recommended in the guidelines such that we have a larger set of peaks to begin with for each replicate. **Peaks were then sorted and only the top 100,000 peaks are kept**. _You do NOT NEED TO RUN this code, we have already generated narrowPeak files for you!_

```
###DO NOT RUN THIS CODE###

# Call peaks using liberal cutoffs
macs2 callpeak -t treatFile.bam -c inputFile.bam -f BAM -g hs -n macsDir/NAME_FOR_OUPUT -B -p 1e-3  2> macsDir/NAME_FOR_OUTPUT_macs2.log

#Sort peak by -log10(p-value)
sort -k8,8nr macsDir/NAME_FOR_OUPUT_peaks.narrowPeak | head -n 100000 > macsDir/NAME_FOR_OUPUT_sorted.narrowPeak

```

> Peak callers tested with IDR:
> 
> * SPP - Works out of the box
> * MACS1.4 - DO NOT use with IDR
> * MACS2 - Works well with IDR with occasional problems of too many ties in ranks for low quality ChIP-seq data.
> * HOMER - developers have a detailed pipeline and code (in beta) for IDR analysis with HOMER at https://github.com/karmel/homer-idr 
> * PeakSeq - Run with modified PeakSeq parameters to obtain large number of peaks
> * HotSpot, MOSAiCS, GPS/GEM, …


### Setting up 

The first thing we need to do is load the module to run IDR:

	$  module load seq/idr/2.0.2


Now let's copy over the narrowPeak filesfor each replicate for Nanog and Pou5f1:

	$ cp



