#' Identify potential mobile sRNA molecules
#'
#' @description A function to identify potential mobile sRNA molecules
#' detected in a tissue which originated from a different genotype, usually 
#' related to a chimeric system.
#'
#' @details
#' The function identifies candidate mobile sRNAs, by selecting sRNA clusters
#' mapped to the foreign genome. It does so by either keeping or removing sRNA 
#' mapped to a given genome. To do so, it requires a common pre-fix across 
#' chromosomes of the same genome. See [mobileRNA::RNAmergeGenomes()] for more
#' information.
#'
#' A greater confidence in the mobile sRNA candidates can be achieved by setting 
#' a threshold that considers the number of replicates which contributed to 
#' defining the consensus dicercall (ie. consensus sRNA classification). This
#' parameter filters based on the `DicerCounts` column introduced by the 
#' [mobileRNA::RNAdicercall()] function.  
#' 
#' 
#' **Statistical Analysis**
#' The function also allows for filtering using statistical inference generated
#' from the differential analysis of the total dataset using the function
#' [mobileRNA::RNAanalysis()]. When `statistical=TRUE`, the feature is enabled
#' and selected mobile molecules that meet the adjusted p-value cutoff defined 
#' by `alpha`. 
#'
#' @param input character; must be either "sRNA" or "mRNA" to represent the type
#' of data, required when setting threshold. 
#'  
#' @param data Numeric data frame
#'
#' @param controls Character vector; containing names of control samples.
#'
#' @param genome.ID a character string related to the chromosomes in a 
#' particular genome. A distinguishing feature of the genome of interest or 
#' non-interest in the chromosome name (`chr` column).
#'
#' @param task an option to keep or remove the chromosomes containing the
#' identifying string. To keep the chromosomes with the ID, set task=keep.
#' To remove, set `task="remove"`. As default, task is set to `keep`.
#'
#'
#' @param statistical If TRUE, will undertake statistical filtering based on the
#' a p-value or adjusted p-value threshold stated by `alpha`.Default
#' set at FALSE. Requires presence of columns containing statistical data.
#' In order to filter by the adjusted p-value, a column named `padjusted` must
#' be present. See [mobileRNA::RNAanalysis()] to calculate statistical values.
#'
#' @param alpha numeric; adjusted p-value cutoff as the target FDR for 
#' independent filtering. Default is 0.1. Only mobile molecules with adjusted 
#' p-values equal or lower than specified are returned.
#'
#'@param threshold numeric; set a threshold level. For sRNAseq, this represents
#' filtering by the minimum number of replicates that defined the 
#' consensus dicercall which is stored in the `DicerCounts` column. 
#'
#' @return A data-frame containing candidate mobile sRNAs, which could be 
#' further filtered based on statistical significance and the ability to 
#' by-pass the threshold which determines the number of replicates that 
#' defined the consensus dicercall. 
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
#'                     genome.ID = "B_",
#'                     task = "keep",
#'                     statistical = FALSE)
#'
#'
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom dplyr "filter"
#' @importFrom dplyr "select"
#' @importFrom tidyselect "starts_with"
#' @importFrom dplyr "case_when"
RNAmobile <- function(data, controls, genome.ID,
                      task = NULL,
                      statistical = FALSE,
                      alpha = 0.1, threshold = NULL){
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
      res <- res %>% filter(padjusted <= alpha)
  } 
  if(!is.null(threshold)){
      res <- res %>% filter(!DicerCounts < threshold)
  }  
  return(res)
}
