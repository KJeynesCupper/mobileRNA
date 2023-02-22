#' Differential Expression (DE) Analysis using `DESeq2` or `edgeR`
#'
#' @description This function allows you to compute the differential expression
#' of sRNA clusters. The function allows the choice between analysis with
#' `"DESeq2"` or `"edgeR"`.
#'
#' @param data numeric data frame produced by `RNAimport()` and/or `RNAsubset()`.
#'
#' @param group Vector of the condition (ie. treatment or control) for each
#' sample. Must be stated in the same order as the samples in the `data` file
#' from left to right.
#'
#' @param method The method to calculate normalisation based on
#' the `DESeq2` or `edgeR` package.
#' Must be stated as either `"DESeq2"` or `"edgeR"`.
#'
#' @param dispersionValue numeric value; manual setting of dispersion value
#' which is recommended for analysis in experiments without biological
#' replicates. This is recommend by `edgeR`.
#' @return The differential expression analysis results.
#'
#' @details The analysis allows the users to choose the method which best suits
#' their data. Notably, `DESeq2` cannot compute the analysis when there only
#' one replicate per condition, but, `edgeR` can. Simply set a suitable
#' dispersion value, based on similar data, to use this feature. The dispersion
#' value is the other wise known as the common Biological  squared coefficient
#' of variation. A typical dispersion value is 0.4 for human data sets, 0.1  for
#' data on genetically identical model organisms or 0.01 for technical replicate.
#' See the User’s Guide for the  ‘EdgeR’ package for more details.
#'
#' @examples
#'# sample conditions.
#' groups <- c("Tomato/Eggplant", "Tomato/Eggplant", "Tomato/Eggplant",
#'           "Tomato/Tomato", "Tomato/Tomato", "Tomato/Tomato")
#'
#'
#'## Differential analysis: DEseq2 method
#'# 24-nt sRNA data-set
#'sRNA_24_DESeq2 <- RNAanalysis(data = sRNA_24,
#'                              group = groups,
#'                              method = "DESeq2" )
#'# 2122-nt sRNA data-set
#'sRNA_2122_DESeq2 <- RNAanalysis(data = sRNA_2122,
#'                                group = groups,
#'                                method = "DESeq2")
#'
#'
#'## Differential analysis: edgeR method
#'sRNA_24_edgeR <- RNAanalysis(data = sRNA_24,
#'                             group = groups,
#'                             method = "edgeR" )
#'
#'sRNA_2122_edgeR <- RNAanalysis(data = sRNA_2122 ,
#'                               group = groups,
#'                               method = "edgeR" )
#' @export
#' @importFrom DESeq2 "results"
#' @importFrom DESeq2 "DESeq"
#' @importFrom stats "relevel"
#' @importFrom edgeR "exactTest"
#' @importFrom stats "p.adjust"
#' @importFrom dplyr "select"
#' @importFrom tidyselect "starts_with"
#' @importFrom dplyr "mutate_at"
#' @importFrom tidyr "replace_na"

RNAanalysis <- function(data, group, method = c("edgeR", "DESeq2"),
                        dispersionValue = NULL){
  if (base::missing(data) || !base::inherits(data,  "data.frame")) {
    stop(paste("Please specify a data frame"))
  }
  counts <- data %>% dplyr::select(tidyselect::starts_with("Count"))

  if (base::missing(method) || !method %in% c("edgeR", "DESeq2")) {
    stop(paste("Please specify analysis method", "(\"edgeR\", or \"DESeq2\")"))
  }
  if ( is.null(group) || !base::inherits(group, "character")) {
    stop("group must be an vector of characters to specify the treatment
         and control conditions for each replicate")
  }
  method <- base::match.arg(method)
  if(method == "edgeR"){
    res <-  .edgeR_normalise(counts, group)
    if(is.null(dispersionValue)){
      mean <- rowMeans(res$counts)
      groupConditions <- unique(group)
      comp <- edgeR::exactTest(res, pair= groupConditions)
      comp <- base::cbind(comp$table, padj=stats::p.adjust(comp$table$PValue,
                                                           method="BH"))
      comp$CountMean <- mean
    } else {
      comp <- edgeR::exactTest(res, dispersion=dispersionValue^2)
      comp <- base::cbind(comp$table, padj=stats::p.adjust(comp$table$PValue,
                                                           method="BH"))
    }
    data$CountMean <- comp$CountMean
    data$log2FoldChange <- comp$logFC
    data$pvalue <- comp$PValue
    data$padjusted <- comp$padj
    data$logCPM <- comp$logCPM
  } else
    if (method == "DESeq2"){
      res <- .DESeq_normalise(counts, group)
      comp <- DESeq2::DESeq(res)
      comp <- DESeq2::results(comp)
      data$baseMean <- comp[,"baseMean"]
      data$log2FoldChange <- comp[,"log2FoldChange"]
      data$pvalue <- comp[,"pvalue"]
      data$padjusted <- comp[,"padj"]
    }
  res.df <- data %>%
    dplyr::mutate_at(c('log2FoldChange','padjusted', 'pvalue'),
                     ~tidyr::replace_na(.,0))
  return(res.df)
}


