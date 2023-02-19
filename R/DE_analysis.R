#' Differential Expression (DE) Analysis using DESeq2 or edgeR
#'
#' @description This function allows you to compute the differential expression
#' of siRNA clusters. The function allows the choice between analysis with
#' DESeq2 or edgeR. This should be chosen based on the normalisation method
#' (See function [DE_prepare()]).
#'
#' @param data Columns containing the raw count data for each sample within
#' the dataframe object.
#' @param groupPair A vector of the two conditions in the comparison
#' @param method The method to calculate normalisation based on
#' the `DESeq2` or `edgeR` package.
#' Must be stated as either `"DESeq2"` or `"edgeR"`.
#' @return The differential expression analysis results.
#' @examples
#'## Differential analysis: DESeq2 example (24-nt & 21/22-nt sRNA)
#' require(DESeq2)
#'
#'
#' # 24-nt sRNA data
#' data("sRNA_24_prep_DESeq2")
#' sRNA_24_DE_DESeq2 <- DE_analysis(data = sRNA_24_prep_DESeq2,
#'                            groupPair = c("Tomato/Eggplant","Tomato/Tomato"),
#'                            method = "DESeq2" )
#'# 21/22-nt sRNA data
#'data("sRNA_2122_prep_DESeq2")
#' sRNA_2122_DE_DESeq2 <- DE_analysis(data = sRNA_2122_prep_DESeq2,
#'                            groupPair = c("Tomato/Eggplant","Tomato/Tomato"),
#'                            method = "DESeq2" )
#'
#'
#'
#'## Differential analysis edgeR example (24-nt & 21/22-nt sRNA)
#' require(edgeR)

#' # 24-nt sRNA data
#' data("sRNA_24_prep_edgeR")
#' sRNA_24_DE_edgeR <- DE_analysis(data = sRNA_24_prep_edgeR,
#'                    groupPair = c("Tomato/Eggplant","Tomato/Tomato"),
#'                      method = "edgeR" )
#' # 21/22-nt sRNA data
#' data("sRNA_2122_prep_edgeR")
#' sRNA_2122_DE_edgeR <- DE_analysis(data = sRNA_2122_prep_edgeR,
#'                            groupPair = c("Tomato/Eggplant","Tomato/Tomato"),
#'                            method = "edgeR" )
#'
#'
#'
#'
#'
#' @export
#' @importFrom DESeq2 "results"
#' @importFrom DESeq2 "DESeq"
#' @importFrom stats "relevel"
#' @importFrom edgeR "exactTest"
#' @importFrom stats "p.adjust"

DE_analysis <- function(data, groupPair, method = c("edgeR", "DESeq2")){
  if ( base::missing(data) || !base::inherits(data,
                                              c("DESeqDataSet", "DGEList"))) {
    stop("data must be an object of class DESeqDataSet
         (if using method == DESeq2) or DGEList (if using method == edgeR)")
  }
  if (base::missing(groupPair) || !base::inherits(groupPair, "character")) {
    stop("groupPair must be an vector of characters to specify the treatment
         and control conditions")
  }
  if (base::missing(method) || !method %in% c("edgeR", "DESeq2")) {
    stop(paste("Please specify analysis method", "(\"edgeR\", or \"DESeq2\")"))
  }
  method <- base::match.arg(method)
  if (method == "edgeR"){
    mean <- rowMeans(data$counts)
    comp <- edgeR::exactTest(data, pair= groupPair)
    comp <- base::cbind(comp$table, padj=stats::p.adjust(comp$table$PValue,
                                                         method="BH"))
    comp$CountMean <- mean
  } else
  if (method == "DESeq2"){
    comp <- DESeq2::DESeq(data)
    comp <- DESeq2::results(comp)
  }
  return(comp)
}




