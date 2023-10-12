# mobileRNA 0.99.9

* updates inline with bioconductor requirements 
* Addition of new function called sRNAmapper
* Deletion of RNAloci function
* Amended RNAmergeGenomes and RNAmergeAnnotations
* Improved vignette 
* Removed lazyloading of package data
* RNAanalysis changed name to RNAdifferentialAnalysis

# mobileRNA 0.99.9
## Previous changes 
* Added a `NEWS.md` file to track changes to the package.
* RNAconsensus changed to RNAdicercall
* Alterations to RNAdicercall algorithm, including tie options, altered default tidy method. 
* RNAdicercall introduced new column "DicerCount" and improved the functionality, specifically the `exclude` parameter.  
* RNAmobile introduced new parameter, "threshold"
* Improved RNAsequences selection algorithm to consider a threshold value, and handling ties. 
* Improved error calling on functions
* Added RNAdf2e function
* Improvements to plotSamplePCA and plotHeatmap
* Amended RNAmergeAnnotations function to meet requirement
* Removed unnecessary man files for GFF and FASTA files on remote repo
* Fixed bug in RNAdistribution plot, when sample specific
* Broadened use of RNAattributes function. 
* Updated vignette
* Updated RNAdicercall to allow any dicer-classification (not constricted to 20-24)
* Converted cat() to message() for user information from functions 
* Amended example data 
* Amended CITATION and NEWs file
* Altered examples in RNAmergeAnnotations/Genomes functions to prevent examples saving into users directory


