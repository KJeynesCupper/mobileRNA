#' Define the sRNA consensus for each dicer-derived cluster
#'
#' @description Using the data, the function uses the supplied data frame and
#' adds an additional column stating the consensus sRNA class/type for each
#' dicer-derived cluster.
#'
#'
#' @details
#' The function calculates the consensus sRNA class based on the conditions
#' supplied. Depending on your reasons for analysis, different conditions should
#' be supplied. For instance, if you wish to identify mobile sRNAs in a
#' heterograft condition, where you compare replicates in either a heterograft
#' or control condition, you should supply the names of the replicates in the
#' heterograft condition. This means that the function will draw a consensus of
#' the sRNA class based only on these replicates. This method is suggested
#' because any mobile sRNA from the donor will not be found in the control
#' samples, and hence, any class determination in the control samples for
#' this sRNA cluster would be irrelevant.
#'
#' @param data a data frame object containing sample data where rows
#' represent sRNA dicer-derived clusters, and where columns represent sample
#' data. See [mobileRNA::RNAimport()] to load data, extract the required
#' information for each sample and organise it as required.
#'
#' @param conditions a character vector containing names sample replicates to
#' base the consensus on. Each string should represent a sample name already
#' utilised in the analysis and present in the data frame supplied to the `data`
#' argument.
#'
#' @param tidy use of this argument will remove sRNA clusters with a unknown or
#'  unclassified consensus result. By default, the function will tidy the data
#'  and remove unclassified clusters to remove excess noise in the data set. To
#'  retain background background noise, set \code{tidy=FALSE}.
#'
#'@return A data frame containing all existing columns in the input data object,
#'plus, an additional column labeled \code{sRNA_Consensus} stating the consensus
#'small RNA type/class between 20-24 nucleotides in length.
#'
#'
#'
#' @examples
#'
#'  data("sRNA_data")
#'
#' # define consensus sRNA classes.
#' conditions <- c("TomEgg_1", "TomEgg_2", "TomEgg_3")
#'
#' # Run function to define sRNA class for each cluster.
#' sRNA_data_summary <- RNAconsensus(data = sRNA_data,
#'                                      conditions = conditions,
#'                                     tidy=TRUE)
#'
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom dplyr "mutate"
#' @importFrom dplyr "select"
#' @importFrom dplyr "filter"
#' @importFrom stringr "str_detect"
#' @importFrom Repitools "annoGR2DF"

RNAconsensus <- function(data, conditions, tidy=TRUE) {
  if (base::missing(data)) {
    stop("data is missing. data must be an object of class matrix,
         data.frame, DataFrame")
  }
  if (!base::inherits(data, c("matrix", "data.frame", "DataFrame"))) {
    stop("data must be an object of class matrix, data.frame, DataFrame")
  }
  if (base::missing(conditions) || !base::inherits(conditions, "character")) {
    stop("conditions must be an vector of characters to define conditions to
         define the consensus from.")
  }
  data[is.na(data)] <- "N"
  class_colnames <- c()
  for (i in colnames(data)){
    if (stringr::str_detect(i, "DicerCall_" )){
      class_colnames <- c(class_colnames, i)
    }
  }
  onlyconditions <- base::unique(grep(paste(conditions,collapse="|"),
                             class_colnames, value=TRUE))
  data <- data %>%
    dplyr::mutate(nt_20 = base::rowSums(.[onlyconditions] == "20"))%>%
    dplyr::mutate(nt_21 = base::rowSums(.[onlyconditions] == "21"))%>%
    dplyr::mutate(nt_22 = base::rowSums(.[onlyconditions] == "22"))%>%
    dplyr::mutate(nt_23 = base::rowSums(.[onlyconditions] == "23"))%>%
    dplyr::mutate(nt_24 = base::rowSums(.[onlyconditions] == "24"))%>%
    dplyr::mutate(nt_N = base::rowSums(.[onlyconditions] == "N"))
  t <-grep('^nt', base::names(data))
  new_df <- data %>%
    dplyr::mutate(sRNA_Consensus = base::names(data)[t]
    [max.col(data[t], ties.method = 'first')*NA^(base::rowSums(data[t])==0)])
  new_df <- new_df %>% dplyr::select(-nt_20,-nt_21,-nt_22,-nt_23, -nt_24,-nt_N)

  if(tidy==TRUE){
    new_df_tidy <- new_df %>% dplyr::filter(sRNA_Consensus != "nt_N")
    return(new_df_tidy)
  }
  else{
    return(new_df) }
}



