#' Define the sRNA dicer consensus for each dicer-derived sRNA cluster
#'
#' @description Using the data, the function uses the supplied dataframe and
#' adds an additional column stating the consensus sRNA class for each
#' dicer-derived cluster.
#'
#' @details
#' The function calculates the consensus sRNA class based on the conditions
#' supplied. 
#' 
#' When ties.method = "random", as per default, ties are broken at random. 
#' In this case, the determination of a tie assumes that the entries are 
#' probabilities: there is a relative tolerance of 1e-5, relative to the 
#' largest (in magnitude, omitting infinity) entry in the row.
#' 
#' When ties.method = "exclude", ties between sRNA classification are ruled as 
#' unclassified ("N"). 
#' 
#' To remove excess data noise, `tidy=TRUE` can be used to removed unclassified 
#' ("N") sRNA clsuters, resulting in a reduced dataset size. 
#' 
#' When working with a chimeric system, for example interspecific grafting, 
#' mapping errors can easily be recognised and eliminated. Here, these can be 
#' eliminated by supplying some extra parameter information. State 
#' `chimeric=TRUE` and supply the chromosome identifier of the foreign genome 
#' (ie. not the tissue sample genotype, but the genotype from which any 
#' potential mobile molecules could be traveling from) to the `genome.ID` 
#' parameter & the control condition samples names to the `controls` parameter.  
#' 
#'
#' @param data a data frame object containing sample data where rows
#' represent sRNA dicer-derived clusters, and where columns represent sample
#' data. See [mobileRNA::RNAimport()] to load data, extract the required
#' information for each sample and organise it as required.
#'
#' @param conditions character vector; containing names of sample replicates.
#' Named replicates will be used to calculate dicercall consensus. 
#'Each string should represent a sample name present in the dataframe supplied 
#'to `data` argument.
#'
#' @param tidy use of this argument will remove sRNA clusters with a unknown or
#'  unclassified consensus result. By default, the function will not tidy the 
#'  data. When \code{tidy=TRUE}, unclassified clusters will be to remove from 
#'  the output dataframe in an effort to remove excess background noise. 
#'  
#' @param ties.method a character string specifying how ties are handled, 
#' "exclude" by default. Options include "random", or 
#' "exclude"; see ‘Details’
#' 
#' 
#'@param chimeric logical; state whether system is chimeric: contains multiple 
#'genomes/genotypes. 
#'
#'@param controls character; vector of control condition sample names. 
#'
#'@param genome.ID character; chromosome identifier of foreign genome in chimeric 
#'system
#'
#'@return A data frame containing all existing columns in the input data object,
#'plus, two additional columns of data: 
#'
#'The first column, `DicerCounts` contains the number of replicates which had 
#'a defined dicer-derived sRNA class (based on the `conditions`). This can be 
#'utilised within the \code{RNAmobile} function as a threshold parameter. 
#'
#'The second, labeled `DicerConsensus` states the consensus sRNA class between 
#'20-24 nucleotides in length or "N" if unclassified. 
#'
#'
#' @examples
#'
#'  data("sRNA_data")
#'
#' # define consensus sRNA classes.
#' conditions <- c("heterograft_1", "heterograft_2", "heterograft_3")
#'
#' # Run function to define sRNA class for each cluster.
#' sRNA_data_consensus <- RNAdicercall(data = sRNA_data,
#'                                   conditions = conditions,
#'                                   tidy=TRUE)
#'
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom dplyr "mutate"
#' @importFrom dplyr "select"
#' @importFrom dplyr "filter"
#' @importFrom stringr "str_detect"
#' @importFrom tidyr "replace_na"
RNAdicercall <- function(data, conditions = NULL, ties.method = NULL, 
                         tidy = FALSE, chimeric = FALSE, controls = NULL, 
                         genome.ID = NULL) {
  if (base::missing(data)) {
    stop("data is missing. data must be an object of class matrix, data.frame, 
         DataFrame")
  }
  if (!base::inherits(data, c("matrix", "data.frame", "DataFrame"))) {
    stop("data must be an object of class matrix, data.frame, DataFrame")
  }
  if(is.null(ties.method)){
    ties.method <-  "exclude"
  }
  data[is.na(data)] <- "N"
  
  # remove mapping errors:
  if(chimeric){
    data <- .remove_mapping_errors_V2(data = data,controls = controls, 
                                      genome.ID = genome.ID)
  }
  
  class_colnames <- c()
  for (i in colnames(data)) {
    if (stringr::str_detect(i, "DicerCall_")) {
      class_colnames <- c(class_colnames, i)
    }
  }
  if (!is.null(conditions)) {
    message("Calculating consensus dicercall based on information from select replicates")
    onlyconditions <- base::unique(grep(paste(conditions, collapse = "|"), 
                                        class_colnames, value = TRUE))
  }
  else if (is.null(conditions)) {
    message("Calculating consensus dicercall based on information from all replicates")
    onlyconditions <- class_colnames
  }
  
  other_exclude <- c("20", "21", "22", "23", "24", "N", "NA")
  data <- data %>% 
    dplyr::mutate(nt_20 = rowSums(.[onlyconditions] =="20")) %>% 
    dplyr::mutate(nt_21 = rowSums(.[onlyconditions] == "21")) %>% 
    dplyr::mutate(nt_22 = rowSums(.[onlyconditions] ==  "22")) %>% 
    dplyr::mutate(nt_23 = rowSums(.[onlyconditions] == "23")) %>% 
    dplyr::mutate(nt_24 = rowSums(.[onlyconditions] == "24"))%>% 
    dplyr::mutate(other = rowSums(!sapply(dplyr::select(.,onlyconditions), 
                                          `%in%`, other_exclude)))
  
  # search columns based on location 
  col_q <- grep("^nt", base::names(data))
  col_qp <- grep("^other", base::names(data))
  t <-c(col_q,col_qp)
  
  if (ties.method == "random"){
    message("The consensus dicercall will be choose at random in the case of a tie")
    new_df <- data %>% 
      dplyr::mutate(DicerConsensus = base::names(data)[t]
                    [max.col(data[t], ties.method = "random")* NA^(
                      rowSums(data[t]) ==0)]) %>%
      dplyr::mutate(DicerConsensus = tidyr::replace_na(DicerConsensus, "N")) 
  } else 
    if(ties.method == "exclude"){
      message("The consensus dicercall will be excluded in the case of a tie") 
      new_df <- data
      
      # Initialize result vector
      result <- vector("character", nrow(new_df))
      result[] <- "N"
      dicer_counts <- vector("character", nrow(new_df))
      dicer_counts[] <- "N"
      # For loop to check for two matching non-zero numbers within the same row
      for (i in 1:nrow(new_df)) {
        row_values <- new_df[t][i, ]
        if(rowSums(row_values) == 0){
          classification <- "N"
          dicer_counts_val <- length(onlyconditions)
        } else {
          non_zero_values <- as.numeric(row_values[row_values != 0])
          values_table <- table(non_zero_values)
          max_value <- max(non_zero_values)
          
          # Check if the maximum value is duplicated
          if (!is.na(values_table[max_value]) && values_table[max_value] > 1) {
            dicer_counts_val <- 0
            classification <- "N"
          } else {
            
            classification <- names(row_values)[max.col(row_values)*NA^(
              rowSums(row_values) == 0)]
            count_max <- max.col(row_values)
            dicer_counts_val <- as.numeric(row_values[,count_max])
            if (is.na(classification)){
              classification <- "N"
              dicer_counts_val <- 0
            }
          }
        }
        result[i] <- classification
        dicer_counts[i] <- dicer_counts_val
      }
      new_df$DicerCounts <- as.numeric(dicer_counts)
      new_df$DicerConsensus <- result
    }
  
  # remove calulation columns 
  new_df <- new_df %>% dplyr::select(-nt_20, -nt_21, -nt_22, 
                                     -nt_23, -nt_24, -other)
  # remove nt from output values
  new_df$DicerConsensus <- gsub("^nt_", "", new_df$DicerConsensus)
  
  
  if (tidy) {
    cat("\n")
    message("Removing sRNA clusters with no consensus dicercall...")
    new_df_tidy <- new_df %>% dplyr::filter(DicerConsensus != "N")
    return(new_df_tidy)
  }
    return(new_df)
}
