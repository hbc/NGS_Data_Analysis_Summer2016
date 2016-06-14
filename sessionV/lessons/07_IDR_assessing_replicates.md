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

<img src=../img/idr_samples.png width=500> 


Common methods for handling replicates includes taking overlapping peak calls across replicates and then assessing differences in binding regions. However, these are simple methods that do not employ any statistical testing and so we know little about how robust these peaks truly are.


## Irreproducibility Discovery Rate (IDR)

[IDR](https://sites.google.com/site/anshulkundaje/projects/idr) is a framework developed by Qunhua Li and Peter Bickel's group that **compares a pair of ranked lists of regions/peaks and assigns values that reflect its reproducibility.**

<img src=../img/idr_figure.png width=500> 

It is extensively used by the ENCODE and modENCODE projects and is part of their [ChIP-seq guidelines and standards](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/). It has been established for submission of ChIP-seq data sets and have been constructed based on the historical experiences of ENCODE ChIP-seq data production groups.

### Why IDR?

* IDR avoids choices of initial cutoffs, which are not comparable for different callers 
* IDR does not depend on arbitrary thresholds and so all regions/peaks are considered. 
* It is based on ranks, so does not require the input signals to be calibrated or with a specific fixed scale (only order matters).


### Components of IDR

1) A **correspondence curve**: a graphical representation of matched peaks as you go down the ranked list. Qualitative, not adequate for selecting signals.

<img src=../img/corr_curve.png width=400> 

2) An **inference procedure**: summarizes the proportion of reproducible and irreproducible signals. Quantitative, using a copula mixture model.

> What proportion of identifications have a poor correspondence, i.e. falling into ”noise”?
> How consistent are the identifications before reaching breakdown?

3) **Irreproducible Discovery Rate (IDR)**: Derive a significance value from the inference procedure (#2) in a fashion similar to FDR, and can be used to control the level of irreproducibility rate when selecting signals.
i.e. 0.05 IDR means that peak has a 5% chance of being an irreproducible discovery



### The IDR pipeline




Full dataset was dowloaded from ENCODE:

* Input: https://www.encodeproject.org/experiments/ENCSR000BHL/
* Nanog Replicates: https://www.encodeproject.org/experiments/ENCSR000BMT/
* Pou5f1 Replicates: https://www.encodeproject.org/experiments/ENCSR000BMU/

Load the module:

	$  module load seq/idr/2.0.2


Copy over the full BEDfiles  Nanog and Pou5f1

	$ cp



