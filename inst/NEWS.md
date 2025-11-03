# mobileRNA 1.7.0

* Added Namespace dependency to DESCRIPTION Imports/Depends entries: ‘Seqinfo’

# mobileRNA 1.5.2

* Updated several functions including RNAfeatures. 

# mobileRNA 1.0.15

* Fixed typo in mapRNA

# mobileRNA 1.0.14

* Fixed inaccuracies in README
* Fixed bug in RNAmergeGenomes


# mobileRNA 1.0.13

* Updates citation 

# mobileRNA 1.0.12

* Updated RNAfeatures and plotRNAfeatures including new parameters to
control plotting and representation of genomic features. 
* Updated citation 

# mobileRNA 1.0.11

* Major updated to RNAfeatures. Instead identifies features independently for 
the genome annotation rather than using preset features.  

# mobileRNA 1.0.1

* Version bump following Bioconductor updates.
* Improved the RNAfeatures function due discrepancies. Previously, calculated 
width rather than relative number of features. 
* Improved RNAdistribution function, including fixing a small bug, replaced 
depreciated functioned and improved aesthetic and code structure. 
* Introduced new plotting function plotRNAfeatures. 

# mobileRNA 0.99.23

* Improved RNAfeatures function 

# mobileRNA 0.99.22

* Fixed issues in RNAimport function


# mobileRNA 0.99.20

* Build redo


# mobileRNA 0.99.19

* Missing connective in RNAfeatures when using repeats variable.
* Added additional check for RNAsequences methods. 

# mobileRNA 0.99.17

* Corrected ORCID references for authors 


# mobileRNA 0.99.17

* Amended `RNAmobile()` to ensure removal of non-zero values 
* Improved look of PCA plot 
* Corrected documentation issue and code disparage in `plotHeatmap()`, and improved styling. 
* Included calculation of consensus sequence determination option for `RNAsequence()`


# mobileRNA 0.99.15

* Alteration of vignette and README
* Included clean FASTQ files as example data sets 
* Updated R data objects
* Updated citation & news files 
* Addition of new function called `mapRNA()`
* Deletion of `RNAloci()` and `RNAmean()` functions
* Additional parameters to `RNAimport()` to support mRNA data importation
* Removed parallel computation in `RNAmergeGenomes()`
* Improved documentation of functions and removed inconsistencies. 



# mobileRNA 0.99.14
## Previous changes 
* Added a `NEWS.md` file to track changes to the package.
* `RNAconsensus()` changed to `RNAdicercall()`
* Alterations to `RNAdicercall()` algorithm, including tie options, altered default tidy method. 
* `RNAdicercall()` introduced new column "DicerCount" and improved the functionality, specifically the `exclude` parameter.  
* `RNAmobile()` introduced new parameter, "threshold"
* Improved `RNAsequences()` selection algorithm to consider a threshold value, and handling ties. 
* Improved error calling on functions
* Added `RNAdf2e()` function
* Improvements to `plotSamplePCA()` and `plotHeatmap()`
* Amended `RNAmergeAnnotations` function to meet requirement
* Removed unnecessary man files for GFF and FASTA files on remote repo
* Fixed bug in `RNAdistribution()` plot, when sample specific
* Broadened use of `RNAattributes()` function. 
* Updated vignette
* Updated `RNAdicercall()` to allow any dicer-classification (not constricted to 20-24)
* Converted cat() to message() for user information from functions 
* Amended example data 
* Amended CITATION and NEWs file
* Altered examples in RNAmergeAnnotations/Genomes functions to prevent examples saving into users directory
* updates inline with bioconductor requirements 
* Amended `RNAmergeGenomes()` and `RNAmergeAnnotations()`
* Removed lazy loading of package data
* `RNAanalysis()` changed name to `RNAdifferentialAnalysis()`

