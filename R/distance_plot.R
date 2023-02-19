#' Heatmap to access sample distance as a quality control step
#'
#'A function to draw a simple hierarchical clustered heatmap to observe
#'sample distance.
#' @param data a `DESeqDataSet` or `DGEList` data object depending on
#' the method see [DE_prepare()] for more information.
#' @return A blue scale heatmap illustrating the sample distance
#' @examples
#'
#' ## DESeq2 example: Sample distance using a Heatmap
#'
#' # plot 24-nt sRNA dataset
#' data("sRNA_24_prep_DESeq2")
#' p1 <- distance_plot(sRNA_24_prep_DESeq2)
#'
#' # plot 21/22-nt sRNA dataset
#' data("sRNA_2122_prep_DESeq2")
#' p2 <- distance_plot(sRNA_2122_prep_DESeq2 )
#'
#'
#'
#' \dontrun{
#' ## edgeR example: Sample distance using a Heatmap
#'
#' # plot 24-nt sRNA dataset
#' data("sRNA_24_prep_edgeR")
#' p3 <- distance_plot(sRNA_24_prep_edgeR)
#'
#' # plot 21/22-nt sRNA dataset
#' data("sRNA_2122_prep_edgeR")
#' p4 <- distance_plot(sRNA_2122_prep_edgeR)
#'}
#'
#' @export
#' @importFrom DESeq2 "rlog"
#' @importFrom stats "dist"
#' @importFrom grDevices "colorRamp"
#' @importFrom pheatmap "pheatmap"
#' @importFrom DEFormats "as.DESeqDataSet"
#' @importFrom SummarizedExperiment "assay"
distance_plot <- function(data){
  if (base::missing(data) || !base::inherits(data, c("DESeqDataSet",
                                                     "DGEList"))) {
    stop("data must be an object of class DESeqDataSet, or DGEList.
         See ?distance_plot for more information.")
  }
  message("Checking data")
  if (base::inherits(data, c("DESeqDataSet")) == TRUE){
    message("Calculating distance matrix")
    rld <- DESeq2::rlog(data, blind = TRUE) # log transform the data.
    distance <- stats::dist(t(SummarizedExperiment::assay(rld)))
    distance_matrix <- as.matrix(distance)
    df_sample_names <- names(data$sizeFactor)
    sample_names <- base::sub('Count_', '', df_sample_names)
    rownames(distance_matrix) <- paste(sample_names)
    colnames(distance_matrix) <- NULL
    message("Creating sample distance plot")
    colors <- grDevices::colorRampPalette(
      rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
    plot <-  pheatmap::pheatmap(distance_matrix,
                                clustering_distance_rows = distance,
                                clustering_distance_cols = distance,
                                col = colors)
  } else
    if (base::inherits(data, c("DGEList")) == TRUE){
      message("Calculating distance matrix")
      ob_from_norm <- base::get("edgeR_count_file", envir = .GlobalEnv)
    dgelist_to_dds <- DEFormats::as.DESeqDataSet(ob_from_norm)
    rld <- DESeq2::rlog(dgelist_to_dds, blind = TRUE) # log transform the data.
    df_sample_names <- colnames(data$counts)
    sample_names <- base::sub('Count_', '',df_sample_names)
    distance <- stats::dist(t(SummarizedExperiment::assay(rld)))
    distance_matrix <- as.matrix(distance)
    rownames(distance_matrix) <- paste(sample_names)
    colnames(distance_matrix) <- NULL
    message("Creating sample distance plot")
    colors <- grDevices::colorRampPalette( rev(
      RColorBrewer::brewer.pal(9, "Blues")) )(255)
    plot <-  pheatmap::pheatmap(distance_matrix,
                                clustering_distance_rows = distance,
                                clustering_distance_cols = distance,
                                col = colors)
  }
  return(plot)
}








