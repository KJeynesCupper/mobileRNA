#' PCA for quality control
#'
#' @description A function to draw a PCA plot. The function undertakes rlog
#'  transformation of the data in an unbiased manner (\code{blind=TRUE}).
#'
#' @param data Either a \code{DESeqDataSet} or \code{DGEList} containing
#' organised or normalised counts based on the method of analysis.
#' @return A PCA plot to show sample distance.
#' @examples
#'
#' ## DESeq2 example: Sample distance using a PCA plot
#'
#' # plot 24-nt sRNA dataset
#'  data("sRNA_24_prep_DESeq2")
#' p1 <- PCA_plot(sRNA_24_prep_DESeq2)
#'
#'
#' # plot 21/22-nt sRNA dataset
#'  data("sRNA_2122_prep_DESeq2")
#' p2 <- PCA_plot(sRNA_2122_prep_DESeq2 )
#'
#'
#'
#' \dontrun{
#' ## edgeR example: Sample distance using a PCA plot
#'
#' # plot 24-nt sRNA dataset
#' data("sRNA_24_prep_edgeR")
#' p3 <- PCA_plot(sRNA_24_prep_edgeR)
#'
#' # plot 21/22-nt sRNA dataset
#' data("sRNA_2122_prep_edgeR")
#' p4 <- PCA_plot(sRNA_2122_prep_edgeR)
#'}
#'
#'
#'
#'
#'
#'
#' @export
#' @importFrom DESeq2 "rlog"
#' @importFrom DESeq2 "plotPCA"
#' @importFrom ggplot2 "position_nudge"
#' @importFrom ggrepel "geom_label_repel"
#' @importFrom DEFormats "as.DESeqDataSet"
PCA_plot <- function(data){
  if (base::missing(data) || !base::inherits(data,
                                             c("DESeqDataSet","DGEList"))) {
    stop("data must be an object of class DESeqDataSet, or DGEList.
         See ?PCA_plot for more information.")
  }
  message("Checking data")
  if (base::inherits(data, c("DESeqDataSet")) == TRUE){
  message("Calulating log transformation")
  rld1 <- DESeq2::rlog(data, blind = TRUE)
  pca <- DESeq2::plotPCA(rld1, returnData = TRUE, intgroup = "conditions")
  # use the DEseq plot pca function, store in an object.
  nudge <- ggplot2::position_nudge(y = 1)## change position
  df_sample_names <- base::names(data$sizeFactor)
  sample_names <- base::sub('Count_', '', df_sample_names)
  pca["ID"] <- sample_names # create new column with sample names
  message("Organising principal component analysis")
  X <-DESeq2::plotPCA(rld1, intgroup = "conditions")+
    ggrepel::geom_label_repel(data = pca, ggplot2::aes(label = ID),
                              position = nudge, show.legend = FALSE)
  } else
    if (base::inherits(data, c("DGEList")) == TRUE){
      message("Calulating log transformation")
      ob_from_norm <- base::get("edgeR_count_file", envir = .GlobalEnv)
      dgelist_to_dds <- DEFormats::as.DESeqDataSet(ob_from_norm)
      rld <- DESeq2::rlog(dgelist_to_dds, blind = TRUE) # log transform the data
      df_sample_names <- colnames(data$counts)
      sample_names <- sub('Count_', '',df_sample_names)
      pca <- DESeq2::plotPCA(rld, returnData = TRUE, intgroup = "group")
      # use the DEseq plot pca function, store in an object.
      nudge <- ggplot2::position_nudge(y = 1)## change position
      pca["ID"] <- sample_names # create new column with sample names
      message("Organising principal component analysis")
      X <-DESeq2::plotPCA(rld, intgroup = "group")+
        ggrepel::geom_label_repel(data = pca, ggplot2::aes(label = ID),
                                  position = nudge, show.legend = FALSE)
    }
  return(X)
}



