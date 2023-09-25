#' Extract statistically significant sRNA clusters within a population. 
#'
#' @description sRNA clusters which are more abundant, show a statistically 
#' significant logfc. 
#'
#' @details Filters the data based on statistical significance
#'which has been calculated by the [mobileRNA::RNAdifferentialAnalysis()] 
#'function. This optional features enables the user to select sRNA clusters 
#'which meet a specific p-value or adjusted p-values threshold.
#'
#'Requires either or both columns: `pvalue`, `padjusted` to undertake. 
#'
#'
#' When working with a chimeric system, for example interspecific grafting, 
#' mapping errors can easily be recognised and eliminated. Here, these can be 
#' eliminated by supplying some extra parameter information. State 
#' `chimeric=TRUE` and supply the chromosome identifier of the foreign genome 
#' (ie. not the tissue sample genotype, but the genotype from which any 
#' potential mobile molecules could be traveling from) to the `genome.ID` 
#' parameter & the control condition samples names to the `controls` parameter.  
#' 
#' @param data Numeric data frame
#'
#' @param statistical If TRUE, will undertake statistical filtering based on the
#' a p-value or adjusted p-value threshold stated by `padj` & `p.value`.Default
#' set at FALSE. Requires presence of columns containing statistical data.
#' In order to filter by the adjusted p-value, a column named `padjusted` must
#' be present. Similarly, to filter by the p-value, a column named `pvalue` must
#' be present. See [mobileRNA::RNAdifferentialAnalysis()] to calculate 
#' statistical values.
#'
#' @param padj A user defined numeric value to represent the adjusted p-value
#' threshold to define statistic significance. Defaults set at 0.05.Only mobile
#' molecules with adjusted p-values equal or lower than specified are returned.
#'
#' @param p.value A user defined numeric value to represent the p-value
#' threshold to define statistic significance. There is no default value, set
#' this instead of using an adjusted p-value to filter molecules. Only mobile
#' molecules with p-values equal or lower than specified are returned.
#'
#'@param chimeric logical; state whether system is chimeric: contains multiple 
#'genomes/genotypes. 
#'
#'@param controls character; vector of control condition sample names. 
#'
#'@param genome.ID character; chromosome identifier of foreign genome in chimeric 
#'system
#'
#' @return A refined version of the working dataframe supplied to the function.
#' The function selects sRNA clusters which meet the statistical threshold, 
#' given the statistical analysis has been undertaken using the 
#' [mobileRNA::RNAdifferentialAnalysis()] function.
#' 
#'
#' @examples
#'
#'data("sRNA_data_consensus")
#'
#' ## sample conditions in order within dataframe
#'groups <- c("Heterograft", "Heterograft", "Heterograft",
#'             "Selfgraft", "Selfgraft", "Selfgraft")
#' 
#' 
#' ## Differential analysis using the DEseq2 method 
#'sRNA_DESeq2 <- RNAdifferentialAnalysis(data = sRNA_data_consensus,
#'                           group = groups,
#'                           method = "DESeq2")
#'                           
#'  ## Select significant based on padjusted                        
#' significant_sRNAs <- RNAsignificant(sRNA_DESeq2, chimeric = TRUE, 
#'                                     genome.ID = "SL", 
#'                                     controls = c("selfgraft_1", "selfgraft_2", 
#'                                     "selfgraft_3"))
#'                                  
#'                                  
#'
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom dplyr "filter"

RNAsignificant <- function(data, statistical = FALSE, padj = 0.05,
                           p.value = NULL, chimeric = FALSE, controls = NULL, 
                           genome.ID = NULL){
  if (!base::inherits(data, c("matrix", "data.frame", "DataFrame"))) {
    stop("data must be an object of class matrix, data.frame, DataFrame")
  }
  if(chimeric){
    data <- .remove_mapping_errors_V2(data = data,controls = controls, 
                                      genome.ID = genome.ID)
  }  
    if (is.null(p.value)) {
      res <- data %>% filter(padjusted <= padj)
    } else
      res <- data %>% filter(pvalue <= p.value)
    return(res)
}
