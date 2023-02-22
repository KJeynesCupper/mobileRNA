#' Heatmap to access sample distance as a quality control step
#'
#'A function to draw a simple hierarchical clustered heatmap to observe
#'sample distance.
#' @param data data.frame object containing raw count data
#' @return A blue scale heatmap illustrating the sample distance
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
plotSampleDistance <- function(data){
  message("Checking data")
  data <- as.matrix(data %>% dplyr::select(tidyselect::starts_with("Count")))
  message("Transforming the count data to the log2 scale")
  rld <- DESeq2::rlog(data, blind = TRUE) # log transform the data.
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




