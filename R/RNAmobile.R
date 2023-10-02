#' Identify potential mobile sRNA & mRNA molecules
#'
#' @description A function to identify potential mobile sRNA or mRNA moleculea
#' detected in a tissue which originated from a different genotype, usually 
#' related to a chimeric system.
#'
#' @details
#' **For sRNAseq:**
#' The function identifies candidate mobile sRNAs, by selecting sRNA clusters
#' mapped to the foreign genome. It does so by either keeping or removing sRNA 
#' mapped to a given genome. To do so, it requires a common pre-fix across 
#' chromosomes of the same genome. See [mobileRNA::RNAmergeGenomes()] for more
#' information.
#'
#' A greater confidence in the mobile sRNA candidates can be achieved by setting 
#' a threshold that considers the number of replicates which contributed to 
#' defining the consensus dicer-call (ie. consensus sRNA classification). This
#' parameter filters based on the `DicerCounts` column introduced by the 
#' [mobileRNA::RNAdicercall()] function.  
#' 
#' **For mRNAseq:**
#' The function identifies candidate mobile mRNAs, by selecting mRNA molecules
#' mapped to the foreign genome. It does so by either keeping or removing mRNA 
#' mapped to a given genome. To do so, it requires a common pre-fix across 
#' chromosomes of the same genome. See [mobileRNA::RNAmergeGenomes()] and 
#' [mobileRNA::RNAmergeAnnotations()] for more information and assistance. 
#'
#' A greater confidence in the mobile mRNA candidates can be achieved by setting 
#' a threshold that considers the number of replicates which contained reads 
#' for the mRNA molecule. This parameter filters based on the `SampleCounts` 
#' column introduced by the [mobileRNA::RNAimport()] function. 
#' 
#' **Statistical Analysis**
#' The function also allows for statistical analysis based on the results
#' collect from differential analysis of the total dataset using the function
#' [mobileRNA::RNAanalysis()]. When `statistical=TRUE`, the feature is enabled
#' and selected mobile molecules that meet a specific p-value or 
#' adjusted p-values.
#'
#' @param input character; must be either "sRNA" or "mRNA" to represent the type
#' of data, required when setting theshold. 
#'  
#' @param data Numeric data frame
#'
#' @param controls Character vector; containing names of control samples.
#'
#' @param genome.ID a character string related to the chromosomes in a particular
#' genome. A distinguishing feature of the genome of interest or non-interest in
#' the chromosome name (`chr` column).
#'
#' @param task an option to keep or remove the chromosomes containing the
#' identifying string. To keep the chromosomes with the ID, set task=keep.
#' To remove, set `task="remove"`. As default, task is set to `keep`.
#'
#'
#' @param statistical If TRUE, will undertake statistical filtering based on the
#' a p-value or adjusted p-value threshold stated by `padj` & `p.value`.Default
#' set at FALSE. Requires presence of columns containing statistical data.
#' In order to filter by the adjusted p-value, a column named `padjusted` must
#' be present. Similarly, to filter by the p-value, a column named `pvalue` must
#' be present. See [mobileRNA::RNAanalysis()] to calculate statistical values.
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
#'@param threshold numeric; set a threshold level. For sRNAseq, this represents
#' filtering by the minimum number of replicates that defined the 
#' dicer-consensus which is stored in the `DicerCounts` column. For mRNAseq,
#' this represents the minimum number of replicates which contained reads at the
#' gene locus, this information is storeds in the `SampleCounts` column.
#'
#' @return A data-frame containing candidate mobile sRNAs, which could be 
#' further filtered based on statistical significance and the ability to 
#' by-pass the threshold which determines the number of replicates that 
#' defined the dicer-consensus. 
#'
#' @examples
#'
#'
#'data("sRNA_data_consensus")
#'
#'
#' # vector of control names
#' controls <- c("selfgraft_1", "selfgraft_2" , "selfgraft_3")
#'
#' # Locate potentially mobile sRNA clusters associated to tomato, no
#' # statistical analysis
#' mobile_df1 <- RNAmobile(data =  sRNA_data_consensus,
#'                     controls = controls,
#'                     genome.ID = "SL40",
#'                     task = "keep",
#'                     statistical = FALSE)
#'
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr "filter"
#' @importFrom dplyr "select"
#' @importFrom tidyselect "starts_with"
#' @importFrom dplyr "case_when"
RNAmobile <- function(data, controls, genome.ID,
                      task = NULL,
                      statistical = FALSE,
                      padj = 0.05, threshold = NULL, input = NULL, 
                      p.value = NULL){
  if (!base::inherits(data, c("matrix", "data.frame", "DataFrame"))) {
    stop("data must be an object of class matrix, data.frame, DataFrame")
  }
  if (base::missing(controls) || !base::inherits(controls, "character")) {
    stop(paste("Please specify a character vector storing names of control
               replicates"))
  }
  if (base::missing(genome.ID) || genome.ID %in% "") {
    stop(paste("Please specify a single character string which is present in
               the all the chromosomes within the genome you wish to keep
               or remove"))
  }
  
  x <- data %>%
    dplyr::filter(dplyr::case_when(
      is.null(task) & base::grepl(genome.ID, chr) ~ TRUE,
      task == "remove" & !base::grepl(genome.ID, chr) ~ TRUE,
      task == "keep" & base::grepl(genome.ID, chr) ~ TRUE,
      TRUE ~ FALSE
    ))
  
  res <- .remove_mapping_errors(data = x, controls = controls)
  
  # Remove rows with no counts 
  count_columns <- as.numeric(grep("^Count", names(res)))
  # Identify rows where all values in Count columns are zero
  refined <- res %>% select(all_of(count_columns))
  rows_to_remove <- apply(refined, 1, function(row) all(row == 0))
  
  # Remove rows with all zero values in Count columns
  res <- res[!rows_to_remove, ]
  
  if (statistical) {
    if (is.null(p.value)) {
      res <- res %>% filter(padjusted <= padj)
    } else
      res <- res %>% filter(pvalue <= p.value)
  } 
# thresholds for different 

  if(!is.null(threshold)){
    if (base::missing(input)) {
      stop(paste("Please specify a character vector of either `sRNA` or `mRNA` 
                 to input parameter when using the threshold parameter"))
    }
    if(input == "sRNA"){
    res <- res %>% filter(!DicerCounts < threshold)
  }  
    if(input == "mRNA"){
      res <- res %>% filter(!SampleCounts < threshold)
      }
    }  
  return(res)
}

