#' Sample quality control heatmap plot
#'
#'A function to draw a simple hierarchical clustered heatmap to observe
#'sample distance. Ideal for quality control of sequencing data.
#'
#'
#' @param data data frame; data set containing the raw data produced as
#' output from [RNAlocate::RNAconsensus()] and/or [RNAlocate:RNAsubset()]
#'
#' @param vst Logical; Variance stabilizing transformation. By default, the
#' function uses a regularized log transformation on the data set, however, this
#' will not suit all experimental designs.
#'
#'
#' @return A blue scale heatmap illustrating the sample distance.
#'
#' @details In special conditions, regularized log transformation will not suit
#' the experimental design. For example, an experimental design without
#' replicates. In this instance, it is preferable to change the default setting
#' and switch to a variance stabilizing transformation method (\code{vst=TRUE}).
#'
#'
#' @examples
#'
#' data("sRNA_24")
#' plotSampleDistance(sRNA_24)
#'
#' @export
#' @importFrom DESeq2 "rlog"
#' @importFrom stats "dist"
#' @importFrom grDevices "colorRamp"
#' @importFrom pheatmap "pheatmap"
#' @importFrom dplyr "select"
#' @importFrom RColorBrewer "brewer.pal"
plotSampleDistance <- function(data, vst = FALSE){
  message("Checking data")
  data <- as.matrix(data %>% dplyr::select(tidyselect::starts_with("Count")))

  if(vst == TRUE){
    message("Transforming the count data with a variance stabilizing transformation")
    rld <- DESeq2::varianceStabilizingTransformation(data, blind = TRUE)
    # log transform the data.
  } else
    if(vst == FALSE) {
      message("Transforming the count data to the log2 scale")
      rld <- DESeq2::rlog(data, blind = TRUE) # log transform the data.
    }

  message("Calculating distance matrix")
  sample_names <- colnames(data)
  sample_names <- sub("Count_", "", sample_names)
  distance <- stats::dist(t(rld))
  distance_matrix <- as.matrix(distance)
  rownames(distance_matrix) <- paste(sample_names)
  colnames(distance_matrix) <- NULL
  message("Creating sample distance plot")
  colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9,"Blues")))(255)
  plot <- pheatmap::pheatmap(distance_matrix,
                             clustering_distance_rows = distance,
                             clustering_distance_cols = distance,
                             col = colors)
  return(plot)
}


