RNAlocate 
======================================================================

Overview
--------

RNAlocate is an `R` package that provides a pipeline for the rapid identification of mobile RNA molecules in 
plant graft systems. The tool provides an improved pipeline for pre-processing and analysis of sequencing data. 

Installation
------------

The latest version of the package can be install directly from this GitHub repo:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("KJeynesCupper/RNAlocate")
```

Functions
---------

The `RNAlocate` package offer the following 7 pre-processing functions and  10 analysis functions:

**Pre-processing**
- `createWorkplace()`  
- `chrModify()`
- `chrInfo()`
- `mergeFiles()`
- `checkQuality()`
- `findOverlap()`
- `identifyClusters()`

**Analysis**
- `RNAimport()`
- `RNAconsensus()`
- `RNAsubset()`
- `RNAdistribution()`
- `plotConsensusFrequency()`
- `plotSamplePCA()`
- `plotSampleDistribution()`
- `RNAanalysis()`
- `MobileMolecules()`
- `plotHeatmap()`


Loading test data
-----------------

To simulate the usage of the package, there is a grafting data set, from a 
tomato-eggplant grafting experiment. The raw data and the pre-processed data
can be found here. The raw data, `raw_sRNA_data`, contains a snippet of the full data set, where
the reads are drawn from chromosome 1 in each genome. While, the pre-processed
data set, `sRNA_data`,  contains the full experimental data and has already been processed by 
the `RNAimport()` function. 

To load the data set, use the following command:

``` r
data("<data-object>")
```

...where `"<data-object>"` is one of the previously mentioned data sets.

Getting help
------------

For additional information on each function, please read through the documentation in the `RNAlocate` package by typing the `?` help operator before any of the function names in the package or by using the `help()` function.

For a step-by-step quick start analysis, consider reading the vignette provided with this package:

``` r
vignette("RNAlocate-Quick-Start")
```

For a more in-depth analysis, consider reading the full vignette provided with this package:


``` r
vignette("RNAlocate")
```
------------------------------------------------------------------------

*Last updated:* 22-02-2023
