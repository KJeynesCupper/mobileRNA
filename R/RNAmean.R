#' Calculate mean FPKM across samples
#'
#' @description Using the data, the function calculates the mean FPKM for the
#' selected samples and adds the results to the working data frame.
#'
#' @param data data frame, see [RNAlocate::RNAimport()] to produce an organised
#' object of sample data.
#'
#' @param conditions Character vector; represent sample names.
#'
#'
#'@return A data frame containing all existing columns in the input data object,
#'plus, an additional columns. One containing the FPKM means, and Count means.
#'
#'As default, the FPKM means will be stored in the column named `mean_FPKM`,
#'while the count means are stored in column names `mean_Count`.
#'
#'@importFrom stringr "str_detect"
#'@importFrom dplyr "mutate"
#'@importFrom dplyr "%>%"
#'@export
#' @examples
#' data("sRNA_data")
#' selected_samples <- c("TomEgg_1", "TomEgg_2", "TomEgg_3")
#' means <- RNAmean(data = sRNA_data, conditions = selected_samples)
#'
RNAmean <- function(data, conditions){
  FPKM_colnames <- c()
  for (i in colnames(data)){
    if (stringr::str_detect(i, "FPKM_" )){
      FPKM_colnames <- c(FPKM_colnames, i)
    }
  }

  fpkm <- base::unique(grep(paste(conditions,collapse="|"),
                                      FPKM_colnames, value=TRUE))
count_colnames <- c()
  for (i in colnames(data)){
    if (stringr::str_detect(i, "Count_" )){
      count_colnames <- c(count_colnames, i)
    }
  }

  count <- base::unique(grep(paste(conditions,collapse="|"),
                            count_colnames, value=TRUE))

  data <- data %>%
    dplyr::mutate(mean_FPKM = base::rowMeans(.[fpkm])) %>%
    dplyr::mutate(mean_Count = base::rowMeans(.[count]))
}


