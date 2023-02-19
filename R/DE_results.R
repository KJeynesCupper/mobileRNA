#' Extract Differential Expression (DE) analysis results and add to the summary
#' dataframe.
#'
#' @description A function to add the differential expression result to the
#' dataframe containing all information on sRNA clusters. This function
#' dependents on which method was used to call the DE analysis.
#' @param Summarydata Summary dataframe containing the collective information on the
#'  sRNA clusters.
#' @param analysisResults Object contain the DE analysis results
#' (see [DE_prepare()] and [DE_analysis()] for prior steps in
#' DE analysis)
#' @param method The method used to calculate the DE analysis: either
#' `DESeq2` or  `edgeR`  method.
#' @return Adds extra columns of data to the dataframe. A column containing the
#' Log2 Fold-Change results, called `log2FoldChange`, a column containing
#' the p-values, called `pvalue`, a column containing the p-adjusted value,
#' called `padjusted`, and if using the `edgeR` method, a column containing the
#' log counts per million, called `logCPM`, are added.
#'
#' @export
#' @importFrom dplyr "mutate_at"
#' @importFrom tidyr "replace_na"
#'
#' @examples
#' # load summary data
#' data('sRNA_24')
#' data('sRNA_2122')
#'
#'
#'
#' ## DESeq2 example (24-nt & 21/22-nt sRNA)
#' # 24nt subset
#' data('sRNA_24_DE_DESeq2')
#' res_sRNA_24_DESeq2 <- DE_results(Summarydata = sRNA_24,
#'                                  analysisResults = sRNA_24_DE_DESeq2 ,
#'                                  method = "DESeq2")
#'
#'# 2122nt subset
#'data('sRNA_2122_DE_DESeq2')
#' res_sRNA_2122_DESeq2 <- DE_results(Summarydata = sRNA_2122,
#'                                  analysisResults = sRNA_2122_DE_DESeq2 ,
#'                                  method = "DESeq2")
#'
#'
#' ## edgeR example (24-nt & 21/22-nt sRNA)
#' require(edgeR)
#'
#'# 24nt subset
#' data('sRNA_24_DE_edgeR')
#' res_sRNA_24_edgeR <- DE_results(Summarydata = sRNA_24,
#'                                  analysisResults = sRNA_24_DE_edgeR,
#'                                  method = "edgeR")
#'data('sRNA_2122_DE_edgeR')
#' res_sRNA_2122_edgeR <- DE_results(Summarydata = sRNA_2122,
#'                                  analysisResults = sRNA_2122_DE_edgeR,
#'                                  method = "edgeR")
#'
#'

DE_results <- function(Summarydata,analysisResults, method = c("edgeR", "DESeq2")){
  if (base::missing(Summarydata)) {
    stop("data is missing. data must be an object of class matrix, data.frame,
         DataFrame")
  }
  if (!base::inherits(Summarydata, c("matrix", "data.frame", "DataFrame"))) {
    stop("data must be an object of class matrix, data.frame, DataFrame")
  }
  if (!base::inherits(analysisResults, c("DESeqResults", "data.frame"))) {
    stop("analysisResults must be an object of DESeqResults
         (if using method == DESeq2) or data.frame (if using method == edgeR)")
  }
  if (base::missing(method) || !method %in% c("edgeR", "DESeq2")) {
    stop(paste("Please specify analysis method", "(\"edgeR\", or \"DESeq2\")"))
  }
  if (method == "edgeR"){
    Summarydata$CountMean <- analysisResults$CountMean
    Summarydata$log2FoldChange <- analysisResults$logFC
    Summarydata$pvalue <- analysisResults$PValue
    Summarydata$padjusted <- analysisResults$padj
    Summarydata$logCPM <- analysisResults$logCPM
  } else
  if (method == "DESeq2"){
    Summarydata$baseMean <- analysisResults[,"baseMean"]
    Summarydata$log2FoldChange <- analysisResults[,"log2FoldChange"]
    Summarydata$pvalue <- analysisResults[,"pvalue"]
    Summarydata$padjusted <- analysisResults[,"padj"]
  }
  Summarydata <- Summarydata %>%
    dplyr::mutate_at(c('log2FoldChange','padjusted', 'pvalue'),
                     ~tidyr::replace_na(.,0))
  return(Summarydata)
}



