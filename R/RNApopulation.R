#' Identify unique sRNA populations between treatment and control conditions
#'
#' @description Identify unique sRNA populations between two conditions. 
#'
#' @details
#' The function undertakes two different roles. First, it identified whether 
#' there are unique sRNA clusters within one condition in comparison to another.
#' For instance, if there are additional sRNA clusters in your treatment 
#' which are not found in your control samples. Keep in mind, this function does
#' not consider the genome of origin unlike the [mobileRNA::RNAmobile()] 
#' function. 
#' 
#'The second role involves filtering the data based on statistical significance
#'which has been calculated by the [mobileRNA::RNAanalysis()] function. This 
#'optional features enables the user to select sRNA clusters which meet a 
#'specific p-value or adjusted p-values threshold.
#'
#' @param data Numeric data frame
#'
#' @param conditions Character vector; containing names of samples within 
#' conditions to locate unique sRNA clusters. 
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
#'
#' @return A refined version of the working dataframe supplied to the function.
#' The function selects sRNA clusters which are only found in replicates within 
#' the same condition supplied to the `conditions` argument and not in other 
#' replicates within the dataset. 
#' 
#' Selection of sRNA clusters which meet the statistical threshold can be 
#' included, given the statistical analysis has been undertaken using the 
#' [mobileRNA::RNAanalysis()] function.
#' 
#'
#' @examples
#'
#'data("sRNA_data_consensus")
#'
#' # vector of control names
#' reps <- c("heterograft_1", "heterograft_2" , "heterograft_3")
#' 
#' heterograft_pop <- RNApopulation(data = sRNA_data_consensus, 
#'                                  conditions = reps)
#'                                  
#'                                  
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr "filter"
#' @importFrom dplyr "select"
#' @importFrom tidyselect "starts_with"
#' @importFrom dplyr "case_when"

RNApopulation <- function(data, 
                          conditions,
                          statistical = FALSE,
                          padj = 0.05,
                          p.value = NULL){
  
  if (!base::inherits(data, c("matrix", "data.frame", "DataFrame"))) {
    stop("data must be an object of class matrix, data.frame, DataFrame")
  }
  if (base::missing(conditions) || !base::inherits(conditions, "character")) {
    stop(paste("Please specify a character vector storing names of control
               replicates"))
  }
  
  conditions_cols <-  data %>% dplyr::select(paste0("Count_", conditions))
  opposite_cols <- data %>% 
    dplyr::select(dplyr::starts_with("Count_")) %>%
    dplyr::select(!colnames(conditions_cols))
    
  output <- data[0,]
  for(i in 1:nrow(data)){
    sum_conditions <- sum(stats::na.omit(as.numeric(data[
      i,colnames(conditions_cols)],na.rm=TRUE)))
    sum_opposite <- sum(stats::na.omit(as.numeric(data[
      i,colnames(opposite_cols)],na.rm=TRUE)))
    if(sum_opposite == 0 && sum_conditions > 0){
      output <- rbind(output, data[i,])
      }
  }
    if (statistical) {
      if (is.null(p.value)) {
        res <- output %>% filter(padjusted <= padj)
      } else
        res <- output %>% filter(pvalue <= p.value)
    } else
      if (statistical == FALSE){
        res <- output
      }
    return(res)
  }

