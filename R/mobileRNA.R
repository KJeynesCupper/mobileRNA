#' mobileRNA: Explore RNA mobilome & population-scale changes
#'
#' Uses small RNA or messenger RNA sequencing data in two conditions and 
#' identifies changes in the RNA population. mobileRNA was primarily designed
#' for the identification of a putative RNA mobilome in a chimeric system. For 
#' example, in plant graft systems. As input, mobileRNA takes sRNAseq or mRNAseq
#' fastq files. Output consists of a dataframe with putative differences 
#' between two conditions along with a number of plots. 
#'
#'
#' The most important functions in the \pkg{mobileRNA} are:
#' \describe{
#' \item{\code{\link{RNAmergeGenomes}}}{Merge two genome assembly files (FASTA)}
#'  \item{\code{\link{RNAmergeAnnotations}}}{Merge two genome annotation files
#'  (GFF)}
#'  \item{\code{\link{mapRNA}}}{Pre-processing of sRNAseq and mRNAseq 
#'  (alignment, raw count, cluster analysis)}
#'  \item{\code{\link{RNAimport}}}{Reads the pre-processing report files in 
#'  to a dataframe for all conditions}
#'  \item{\code{\link{RNAdicercall}}}{Calculates the consensus sRNA 
#'  dicercall class}
#'  \item{\code{\link{RNAsubset}}}{Subsets the working data set to select a 
#'  defined sRNA class}
#'  \item{\code{\link{RNAdifferentialAnalysis}}}{Undertakes differential 
#'  analysis using either the edgeR or DESeq2 method for data structure}
#'  \item{\code{\link{RNAmobile}}}{Identifies mobile candidates in chimeric 
#'  systems}
#'  \item{\code{\link{RNApopulation}}}{Selects unique RNAs within the defined
#'  sample replicates}
#'  \item{\code{\link{RNAsummary}}}{Summarise the differential abundance of RNAs}
#'   \item{\code{\link{RNAreorder}}}{Reorder the data frame for differential 
#'   analysis, ensuring control verse treatment comparison}
#'  \item{\code{\link{RNAsequences}}}{Defines the consensus RNA sequence for a 
#'  sRNA cluster}
#'  \item{\code{\link{RNAattributes}}}{Overlaps the sRNA clusters with a GFF
#'  annotation file and adds overlapping features to the sRNA cluster 
#'  information}
#'  \item{\code{\link{RNAdistribution}}}{Plot the distribution of sRNA 
#'  lengths/classes}
#'  \item{\code{\link{plotHeatmap}}}{Plots a heatmap of log2-transformed 
#'  normalised RPM/FPKM values}
#'  \item{\code{\link{plotSampleDistance}}}{Plots a sample distance heatmap for
#'  quality control}
#'  \item{\code{\link{plotSamplePCA}}}{Plots a PCA plot, customise ratio, 
#'  colours and shapes}
#'  \item{\code{\link{RNAfeatures}}}{Overlaps the sRNA clusters with a GFF
#'  annotation file and calculates the percentage of clusters which overlap with 
#'  specific genomic features}
#'  \item{\code{\link{RNAdf2se}}}{Convert `mobileRNA` dataframe to 
#'  a SummarizedExperiment object }
#' }
#'
#' @author
#' Katie Jeynes-Cupper \email{kej031@@student.bham.ac.uk},
#' Marco Catoni \email{m.catoni@@bham.ac.uk}
#' Maintainer: Katie Jeynes-Cupper \email{kej031@@student.bham.ac.uk}
#' @name mobileRNA
#' @docType package
#' @seealso See \code{vignette("mobileRNA", package = "mobieRNA")} for an 
#' overview of the package.
NULL


#' @name sRNA_data
#' @title sRNA_data: simulated data for biological replicates
#' @docType data
#' @usage data(sRNA_data)
#' @description Simulated sRNAseq dataset
#' @return Dataframe in global environment 
#' @details Simulates data is taken from eggplant and tomato sRNAseq samples and
#' created to simulate to movement of sRNA molecules from an Tomato rootstock to
#' an Eggplant Scion. Three Eggplant replicates were spiked with the same 150
#' tomato sRNA clusters, and named "heterograft_" 1 to 3. The analysis compares
#' these heterografts to three Eggplant self-grafts which are the original
#'  un-spiked Eggplant replicates, called "selfgraft_" 1 to 3. 
#'  
#'  This data was imported and organised by the [RNAimport()] function. 
#' @examples
#'  data("sRNA_data")
NULL


#' @name mRNA_data
#' @title mRNA_data: simulated messenger RNA data for biological replicates
#' @docType data
#' @usage data(mRNA_data)
#' @description Simulated mRNAseq dataset
#' @return Dataframe in global environment 
#' @details Simulates data is taken from eggplant and tomato mRNAseq samples and
#' created to simulate to movement of mRNA molecules from an Tomato rootstock to
#' an Eggplant Scion. Two Eggplant replicates were spiked with the same 150
#' tomato mRNA clusters, and named "heterograft_" 1 to 2. The analysis compares
#' these heterografts to two Eggplant self-grafts which are the original
#'  un-spiked Eggplant replicates, called "selfgraft_" 1 to 2. 
#'  
#'  This data was imported and organised by the [RNAimport()] function. 
#' @examples
#'  data("mRNA_data")
NULL

