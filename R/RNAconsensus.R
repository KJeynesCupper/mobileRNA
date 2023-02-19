#' Define the siRNA consensus for each dicer-derived cluster
#'
#' @description Using the data, the function creates a data frame with an
#' additional column stating the consensus siRNA class/type for each
#' dicer-derived cluster `sRNA_Consensus`.
#' @param data data frame or GRanges object containing sample data. See
#' [sample_table()] to produce an organised GRanges object of sample data.
#' @param conditions Vector containing names of the heterografted conditions.
#' Or, samples you wish to draw the consensus class from in the analysis.
#' @param tidy Whether to remove clusters with an unknown or unclassified
#' consensus result. `tidy=TRUE` to remove unclassified cluster, or
#' \code{tidy=FALSE} to not remove background noise.
#' It is preferable to remove the excess noise in data.
#'
#'@return A data frame containing all existing columns in the input data object,
#'plus, an additional column labeled `sRNA_Consensus` stating the consensus
#'small RNA type/class between 20-24 nucleotides in length.
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
  }
  else{
    return(new_df) }
}

