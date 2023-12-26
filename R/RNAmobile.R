#' Identify potential mobile RNA molecules
#'
#' @description A function to identify the putative sRNA or mRNA mobilome. 
#' 
#' 
#' @details
#' The function identifies candidate mobile sRNAs or mRNAs, by selecting those
#' that mapped to the genome of mobile molecules (Genome B). It does so by
#' either keeping or removing those mapped to a given genome. To do so, it 
#' requires a common pre-fix across chromosomes of the same genome. See 
#' [mobileRNA::RNAmergeGenomes()] for more information.
#'
#' **For sRNAseq:**
#' A greater confidence in the mobile sRNA candidates can be achieved by setting 
#' a threshold that considers the number of replicates which contributed to 
#' defining the consensus dicercall (ie. consensus sRNA classification). This
#' parameter filters based on the `DicerCounts` column introduced by the 
#' [mobileRNA::RNAdicercall()] function.  
#' 
#' **For mRNAseq:**
#' A greater confidence in the mobile mRNA candidates can be achieved by setting 
#' a threshold that considers the number of replicates which contained reads 
#' for the mRNA molecule. This parameter filters based on the `SampleCounts` 
#' column introduced by the [mobileRNA::RNAimport()] function. 
#' 
#' **Statistical Analysis**
#' The function also allows for filtering using statistical inference generated
#' from the differential analysis of the total dataset using the function
#' [mobileRNA::RNAdifferentialAnalysis()]. When `statistical=TRUE`, the feature 
#' is enabled and selected mobile molecules that meet the adjusted p-value 
#' cutoff defined by `alpha`. 
#'
#' @param input character; must be either "sRNA" or "mRNA" to represent the type
#' of data, required when setting threshold. 
#'  
#' @param data data.frame; generated through the \pkg{mobileRNA} method. 
#'
#' @param controls character vector; containing names of control samples.
#'
#' @param genome.ID character; string related to the chromosomes in a 
#' particular genome. A distinguishing feature of the genome of interest or 
#' non-interest in the chromosome name (`chr` column).
#'
#' @param task character; string to set the method to keep or remove the
#' chromosomes containing the identifying string. To keep the chromosomes with 
#' the ID, set task=keep. To remove, set `task="remove"`. As default, task is 
#' set to `keep`.
#'
#'
#' @param statistical If TRUE, will undertake statistical filtering based on the
#' a p-value or adjusted p-value threshold stated by `alpha`.Default
#' set at FALSE. Requires presence of columns containing statistical data.
#' In order to filter by the adjusted p-value, a column named `padjusted` must
#' be present. See [mobileRNA::RNAdifferentialAnalysis()] to calculate 
#' statistical values.
#'
#' @param alpha numeric; adjusted p-value cutoff as the target FDR for 
#' independent filtering. Default is 0.1. Only mobile molecules with adjusted 
#' p-values equal or lower than specified are returned.
#'
#'@param threshold numeric; set a threshold level. For sRNA analysis, this 
#'represents filtering by the minimum number of replicates that defined the 
#' consensus dicercall which is stored in the `DicerCounts` column. While, 
#' for mRNA analysis this represents the number of replicates which contained 
#' reads for the mRNA molecule which is stored in the `SampleCounts` column. 
#'
#' @return A data-frame containing candidate mobile sRNAs or mRNAs, which could 
#' have been further filtered based on statistical significance and the ability 
#' to by-pass the thresholds which determine the number of replicates that 
#' defined the consensus dicercall (sRNA) or contributed to reads counts (mRNA). 
#'
#' @examples
#'
#'
#'data("sRNA_data")
#'
#'
#' # vector of control names
#' controls <- c("selfgraft_1", "selfgraft_2" , "selfgraft_3")
#'
#' # Locate potentially mobile sRNA clusters associated to tomato, no
#' # statistical analysis
#' mobile_df1 <- RNAmobile(input = "sRNA", data =  sRNA_data,
#' controls = controls, genome.ID = "B_", task = "keep", statistical = FALSE)
#'
#'
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom dplyr "filter"
#' @importFrom dplyr "select"
#' @importFrom tidyselect "starts_with"
#' @importFrom dplyr "case_when"
RNAmobile <- function(input = c("sRNA", "mRNA"), data, controls, genome.ID,
                      task = NULL, statistical = FALSE, alpha = 0.1, 
                      threshold = NULL){
  if (base::missing(input)) {
    stop(paste("Please specify a character vector of either `sRNA` or `mRNA` 
                 to input parameter."))
  }
  
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
  
  y <- data %>%
    dplyr::filter(dplyr::case_when(
      is.null(task) & base::grepl(genome.ID, chr) ~ TRUE,
      task == "remove" & !base::grepl(genome.ID, chr) ~ TRUE,
      task == "keep" & base::grepl(genome.ID, chr) ~ TRUE,
      TRUE ~ FALSE
    ))
  res <- .remove_mapping_errors(data = y, controls = controls)
  if (statistical) {
      res <- res %>% dplyr::filter(padjusted <= alpha)
  } 
  if(!is.null(threshold)){
    if(input == "sRNA"){
      res <- res %>% filter(!DicerCounts < threshold)
    }  
    if(input == "mRNA"){
      res <- res %>% filter(!SampleCounts < threshold)
    }
  }  
  return(res)
}
