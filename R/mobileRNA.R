#' mobileRNA: Explore candidate mobile sRNAs & sRNAs population-scale changes
#'
#' Uses small RNA sequencing data in two conditions and identifies
#' changes in the small RNA population, but mainly, identify foriegn mobile 
#' small RNAs in chimeric systems. For example, in plant graft systems. The 
#' input is the Results report files produced by ShortStack and the output is a 
#' dataframe of small RNA clusters, including a consensus dicercall and 
#' differential analysis results. Output can also consist of the consensus RNA
#' sequence of small RNA clusters within the determined population. Includes
#' quality control and plotting functions. 
#'
#' The most important functions in the \pkg{mobileRNA} are:
#' \describe{
#' \item{\code{\link{RNAloci}}}{reads the ShortStack report files to identify 
#'  all de novo sRNA cluster loci across the experimental design}
#'  \item{\code{\link{RNAimport}}}{reads the ShortStack Results report files in 
#'  to a dataframe for all conditions}
#'  \item{\code{\link{RNAdicercall}}}{calculates the consensus sRNA 
#'  dicercall class}
#'  \item{\code{\link{RNAsubset}}}{subsets the working data set to select a 
#'  defined sRNA class}
#'  \item{\code{\link{RNAdifferentialAnalysis}}}{undertakes differential 
#'  analysis using either the edgeR or DESeq2 method for data structure}
#'  \item{\code{\link{RNAmobile}}}{identifies candidates to be potentially 
#'  mobile sRNA molecules in a chimeric system}
#'  \item{\code{\link{RNApopulation}}}{selects the sRNAs unique to the defined
#'  sample replicates}
#'  \item{\code{\link{RNAsummary}}}{Summarise the abundance of sRNA clusters}
#'   \item{\code{\link{RNAreorder}}}{Reorder the data frame for differential 
#'   analysis, ensuring control verse treatment comparison}
#'  \item{\code{\link{RNAsequences}}}{defines the consensus RNA sequence for a 
#'  sRNA cluster}
#'  \item{\code{\link{RNAattributes}}}{overlaps the sRNA clusters with a GFF
#'  annotation file and adds overlapping features to the sRNA cluster 
#'  information}
#'  \item{\code{\link{RNAdistribution}}}{Plot the distribution of sRNA 
#'  lengths/classes}
#'  \item{\code{\link{plotHeatmap}}}{plots a heatmap of log2-transformed 
#'  normalised RPM/FPKM values}
#'  \item{\code{\link{plotSampleDistance}}}{plots a sample distance heatmap for
#'  quality control}
#'  \item{\code{\link{plotSamplePCA}}}{plots a PCA plot, customise ratio, 
#'  colours and shapes}
#'  \item{\code{\link{RNAfeatures}}}{overlaps the sRNA clusters with a GFF
#'  annotation file and calculates the percentage of clusters which overlap with 
#'  specific genomic features}
#'  \item{\code{\link{RNAmean}}}{calculate mean RPM and mean raw count across 
#'  given samples}
#'  \item{\code{\link{RNAmergeGenomes}}}{Merge two FASTA genome assembly files}
#'  \item{\code{\link{RNAmergeAnnotations}}}{Merge two GFF genome annotation 
#'  files}
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

#' @name sRNA_data_dicercall
#' @title sRNA_data_dicercall: simulated data for biological replicates
#' @docType data
#' @keywords dataset
#' @usage data(sRNA_data_dicercall)
#' @description Simulated sRNAseq dataset with RNA consensus
#' @details Simulates data is taken from eggplant and tomato sRNAseq samples and
#' created to simulate to movement of sRNA molecules from an Tomato rootstock to
#' an Eggplant Scion. Three Eggplant replicates were spiked with the same 150
#' tomato sRNA clusters, and named "heterograft_" 1 to 3. The analysis compares
#' these heterografts to three Eggplant self-grafts which are the original
#' un-spiked Eggplant replicates,  called "selfgraft_" 1 to 3.
#' 
#' This data is a further development of the `sRNA_data` data set, where the 
#' [RNAdicercall()] function has been employed. 
#' @return Dataframe in global environment 
#' @examples
#'  data("sRNA_data_dicercall")
NULL

#' @name sRNA_data_mobile
#' @title sRNA_data_mobile: simulated data for biological replicates
#' @docType data
#' @keywords dataset
#' @usage data(sRNA_data_mobile)
#' @description Simulated sRNAseq dataset - potentially mobile RNAs
#' @return Dataframe in global environment 
#' @details Simulates data is taken from eggplant and tomato sRNAseq samples and
#' created to simulate to movement of sRNA molecules from an Tomato rootstock to
#' an Eggplant Scion. Three Eggplant replicates were spiked with the same 150
#' tomato sRNA clusters, and named "heterograft_" 1 to 3. The analysis compares
#' these heterografts to three Eggplant self-grafts which are the original
#' un-spiked Eggplant replicates, called "selfgraft_" 1 to 3. 
#' 
#' This data is a further development of the `sRNA_data` data set, where the
#' candidate mobile sRNAs have been selected. 
#' @examples
#'  data("sRNA_data_mobile")
NULL

