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

The `RNAlocate` package offer the following 7 pre-processing functions and  11 analysis functions:

**Pre-processing**
- `dir_relocate()`
- `chr_checker()`
- `chr_info()`
- `merge_files()`
- `quality_check()`
- `find_overlap()`

**Analysis**
- `RNAimport()`
- `RNAconsensus()`
- `RNAsubset()`
- `RNAdistribution()`
- `frequency_plot()`
- `DE_prepare()`
- `PCA_plot()`
- `distance_plot()`
- `DE_analysis()`
- `DE_results()`
- `mobile_molecules()`
- `heatmap_plot()`


Loading test data
-----------------

To simulate the usage of the package, there are two grafting data sets. Both of the data sets
were produced from the same raw data, yet, the only difference is the method of calling differential expression. 
The package has the flexibility of the user to choose either the `DESeq2` or `edgeR` software. 
Therefore, up to this point in the analysis, the data is identical. 
-   `sRNA-DESeq-data`
-  `sRNA-edgeR-data`


To load these data sets, use the following command:

``` r
data("<object-type>")
```

...where `"<object-type>"` is one of the previously mentioned data sets.

Getting help
------------

For additional information on each function, please read through the documentation in the `RNAlocate` package by typing the `?` help operator before any of the function names in the package or by using the `help()` function.

For a step-by-step quick start analysis, consider reading the vignette provided with this package:

``` r
vignette("RNAlocate:Quick _Start")
```

For a more indepth analysis, consider reading the full vignette provided with this package:


``` r
vignette("RNAlocate")
```
------------------------------------------------------------------------

*Last updated:* 02-02-2023
