#' Calculate FPKM across samples
#'
#' @description Using the data, the function calculates the mean RPM and Counts
#' across all or specific samples.
#'
#' @param data data frame, see [RNAlocate::RNAimport()] to produce an organised
#' object of sample data.
#'
#' @param conditions Character vector; represent sample names. If not supplied,
#' function will calculate means across all samples in dataframe.
#'
#'
#'@return A data frame containing all existing columns in the input data object,
#'plus, an additional columns.
#'
#'As default, the FPKM means will be stored in the column named `mean_RPM`,
#'while the count means are stored in column names `mean_Count`.
#'
#'@importFrom stringr "str_detect"
#'@importFrom dplyr "mutate"
#'@importFrom dplyr "select"
#'@importFrom dplyr "%>%"
#'@export
#' @examples
#' data("sRNA_data")
#' # across specific samples
#' selected_samples <- c("TomEgg_1", "TomEgg_2", "TomEgg_3")
#' means <- RNAmean(data = sRNA_data, conditions = selected_samples)
#' # for all samples
#'means <- RNAmean(data = sRNA_data)
RNAmean <- function(data, conditions = NULL){
  RPM_cols <- dplyr::select(data, starts_with("RPM_"))
  count_cols <- dplyr::select(data, starts_with("Count_"))

  if (!is.null(conditions)){
    RPM <- base::unique(grep(paste(conditions,collapse="|"),
                             RPM_cols, value=TRUE))
    count <- base::unique(grep(paste(conditions,collapse="|"),
                               count_cols, value=TRUE))
    data <- data %>%
      dplyr::mutate(mean_RPM = base::rowMeans(.[RPM])) %>%
      dplyr::mutate(mean_Count = base::rowMeans(.[count]))

  }
  else
    if (is.null(conditions)){
      data <- data %>%
        dplyr::mutate(mean_RPM = base::rowMeans(.[RPM_cols])) %>%
        dplyr::mutate(mean_Count = base::rowMeans(.[count_cols]))
    }
  return(data)
}



