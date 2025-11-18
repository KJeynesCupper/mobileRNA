#' mobileRNA: Explore RNA mobilome & population-scale changes
#'
#' Uses small RNA or messenger RNA sequencing data in two conditions and
#' identifies changes in the RNA population. `mobileRNA` was primarily designed
#' for the identification of a putative RNA mobilome in a chimeric system. For
#' example, in plant graft systems. As input, `mobileRNA` takes sRNAseq or
#' mRNAseq fastq files. Output consists of a data frame with putative
#' differences between two conditions along with a number of plots.
#'
#'
#' The most important functions in the \pkg{mobileRNA} are:
#' \describe{
#' \item{\code{\link{RNAmergeGenomes}}}{Merge two genome assembly files (FASTA).}
#'  \item{\code{\link{RNAmergeAnnotations}}}{Merge two genome annotation files
#'  (GFF).}
#'  \item{\code{\link{mapRNA}}}{Pre-processing of sRNAseq and mRNAseq
#'  (alignment, raw count, cluster analysis).}
#'  \item{\code{\link{RNAimport}}}{Reads the pre-processing report files into a
#'  dataframe for all conditions}
#'  \item{\code{\link{RNAdicercall}}}{Calculates the consensus sRNA
#'  dicercall class.}
#'  \item{\code{\link{RNAsubset}}}{Subsets the data set based on the sRNA class.}
#'  \item{\code{\link{RNAdifferentialAnalysis}}}{Undertakes differential
#'  analysis with either the edgeR or DESeq2 method.}
#'  \item{\code{\link{RNAmobile}}}{Identify putative RNA molecules produced by
#'  the non-tissue sample genome}
#'  \item{\code{\link{RNApopulation}}}{Identify gained/lost RNA populations
#'  between treatment and control conditions.}
#'  \item{\code{\link{RNAsummary}}}{Summarise the differential abundance of
#'  RNAs.}
#'   \item{\code{\link{RNAreorder}}}{Reorder the data frame for differential
#'   analysis, ensuring control verse treatment comparison.}
#'  \item{\code{\link{RNAsequences}}}{Extract RNA sequence from sRNA clusters.}
#'
#'  \item{\code{\link{RNAattributes}}}{Overlap the genomic features related to
#'  the sRNA clusters.}
#'  \item{\code{\link{RNAdistribution}}}{Plot the distribution of sRNA classes
#'  based on nucleotide length.}
#'  \item{\code{\link{plotHeatmap}}}{Heatmap of log-transformed normalization
#'  data.}
#'  \item{\code{\link{plotSampleDistance}}}{Plots a sample distance heatmap for
#'  quality control.}
#'  \item{\code{\link{plotSamplePCA}}}{Plots a PCA plot, customize ratio,
#'  colours and shapes.}
#'  \item{\code{\link{RNAfeatures}}}{Summarise the distribution of sRNA clusters
#'  across genomic features.}
#'  \item{\code{\link{RNAdf2se}}}{Convert a `mobileRNA` dataframe to a
#'  SummarizedExperiment object. }
#' }
#'
#' @rdname mobileRNA-package
#' @name mobileRNA
#' @keywords internal
#' @aliases mobileRNA-package mobileRNA
#' @docType package
#' @author
#' Katie Jeynes-Cupper \email{kejc@@illinois.edu},
#' Marco Catoni \email{m.catoni@@bham.ac.uk}
#' Maintainer: Katie Jeynes-Cupper \email{kejc@@illinois.edu}
#' @seealso See \code{vignette("mobileRNA", package = "mobieRNA")} for an
#' overview of the package.
"_PACKAGE"
## usethis namespace: start
## usethis namespace: end
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

