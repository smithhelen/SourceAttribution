Source Attribution
================

-   <a href="#overview" id="toc-overview">Overview</a>
-   <a href="#contents-of-this-repository"
    id="toc-contents-of-this-repository">Contents of this repository</a>

### Overview


This repository is to accompany the thesis **"Source Attribution models using Whole Genome Sequencing data"**, 
by *Helen L Smith*
under the supervisory team of *JC Marshall, NP French, PJ Biggs, and ANH Smith*, Massey University, 2025.

This thesis tackles the issue of absent levels in random forest predictive models in the context of a source attribution study of *Campylobacter jejuni* and *C. coli* in New Zealand.

It details three new methods for encoding categorical variables for random forest predictive models which are unbiased in the presence of absent levels 
- the CA-unbiased-encoding method
- the PCO-encoding method, and
- the CAP-encoding method.

It also adapts the MDA (mean decrease in accuracy, or permutation) method of variable importance to be independent of the out-of-bag data and therefore unbiased when target-based methods of encoding are used.

This repository contains all the code for these methods, as well as a file demonstrating how to apply the methods to the source attribution of *Campylobacter* species, including the calculation of uncertainty intervals.

### Contents of this repository

The key contents are organised as follows:

-   Source Attribution
    -   README.Rmd
    -   SourceAttribution.Rproj
    -   data
        -   SourceAttribution.Rdata
        -   SeqDat.Rdata
    -   files
        -   sample_run.R
    -   methods
        -   ca_unbiased.R
        -   pco.R
        -   cap.R
        -   hamming.R
        -   residualise.R
        -   helpers.R
        -   vi_independent_holdout.R
        -   recipe_ca0.R
        -   recipe_pco.R
        -   recipe_cap.R


The R project is called `SourceAttribution.Rproj`.

The directory `methods` is where all of the functions for this project are stored.

-   `ca_unbiased.R` contains the code for the adapted ca method where
    observations with absent levels are directed according to the *a
    priori* hypothesis of equal class distribution
-   `pco.R` contains the code for our new method where observations with
    absent levels are directed through the random forest according to
    their *similarity* to other levels
-   `cap.R` contains the code for our new method where observations with
    absent levels are directed through the random forest according to
    their *similarity* to other levels under the constraint that the
    direction of similarity must be in the direction of greatest class
    difference
-   `hamming.R` for generating a matrix of Hamming distances between
    levels of a factor
-   `helpers.R` contains code for functions used internally in the other
    methods
-   `recipe_ca0.R` contains the ca unbiased method in tidymodels form
-   `recipe_pco.R` contains the pco method in tidymodels form        
-   `recipe_cap.R` contains the cap method in tidymodels form        

The directory `files` contains one file 
-   `sample_run.R` which includes an example of running the methods and generating source attribution estimates


### Data for exploration

The data to accompany the example file is a small dataset for the
purpose of demonstrating the methods in this respository

There are two files:
1. `SeqDat.R` contains the aligned nucleotide
sequencing data for 10 genes from each of 40 *Camplyobacter* isolates
collected from four sources (Sheep, Cattle, Chicken, and Human)
2. `CapItOff.R` contains the corresponding allelic information
for the aligned sequences.
