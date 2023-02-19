#' Prepare data for differential expression analysis with
#' either DESeq2 or edgeR
#'
#' @description A function utilizing `DESeq2` or `edgeR` to prepare
#' data for differential expression analysis. The two methods of analysis
#' require different inputs for analysis, the main difference being that DESeq2
#' expects un-normalised counts whilst `edgeR` expects normalised counts.
#'
#'
#'
#'
#' @param data Columns containing the raw count data for each sample within the
#' Granges object.
#' @param conditions Vector of the condition for each sample, must be stated
#' in the same order as the samples in
#' @param method The method to calculate normalisation based on the
#' `DESeq2` or `edgeR` package. Must be stated as either `"DESeq2"`
#' or `"edgeR"`
#' @return A `DESeqDataSet` object or `DGEList` containing
#' organised or normalised count data for each sample.
#'
#' When using the `edgeR` method, the function produces an additional
#' variable which is stored in the global environment. This object is called
#' `edgeR_count_file` and contains the `DGEList` object produced by
#' `edgeR` containing the un-normalised counts. This variable is required
#' for functions [PCA_plot()] and [distance_plot()].
#' @examples
#'
#'
#' data(sRNA_24)
#' data(sRNA_2122)
#'
#'
#'  ## DESeq2 example (24-nt & 21/22-nt sRNA)
#'
#' # sample conditions.
#' group <- c("Tomato/Eggplant", "Tomato/Eggplant", "Tomato/Eggplant",
#'          "Tomato/Tomato", "Tomato/Tomato", "Tomato/Tomato")
#'
#' # Normalise the 24-nt dataset
#' sRNA_24_prep_DESeq2 <- DE_prepare(sRNA_24, group, method = "DESeq2" )
#'
#' # Normalise the21/22-nt dataset
#' sRNA_2122_norm_DESeq2 <- DE_prepare(sRNA_2122, group, method = "DESeq2" )
#'
#'
#'
#'
#'
#' ## edgeR example (24-nt & 21/22-nt sRNA)
#'
#' # sample conditions.
#' group <- c("Tomato/Eggplant", "Tomato/Eggplant", "Tomato/Eggplant",
#'          "Tomato/Tomato", "Tomato/Tomato", "Tomato/Tomato")
#'
#' # Normalise the 24-nt dataset
#' sRNA_24_norm_edgeR <- DE_prepare(sRNA_24, group, method = "edgeR" )
#'
#'
#' # Normalise the21/22-nt dataset
#' sRNA_2122_norm_edgeR <- DE_prepare(sRNA_2122, group, method = "edgeR" )

#'
#'
#'
#'
#' @export
#' @importFrom dplyr "select"
#' @importFrom dplyr "%>%"
#' @importFrom tidyselect "starts_with"
DE_prepare <- function(data, conditions, method = c("edgeR", "DESeq2")){
  data <- data %>% dplyr::select(tidyselect::starts_with("Count"))
  if (base::missing(method) || !method %in% c("edgeR", "DESeq2")) {
    stop(paste("Please specify analysis method", "(\"edgeR\", or \"DESeq2\")"))
  }
    method <- base::match.arg(method)
  if(method == "edgeR"){
    res <-  .edgeR_normalise(data, conditions)
  } else
    if (method == "DESeq2"){
      res <- .DESeq_normalise(data, conditions)
  }
  return(res)
}
