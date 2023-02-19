#' Heatmap using hierarchical clustering
#'
#' @description Plot a heatmap with hierarchical clustering via an rlog
#' transformation of FPKM data and euclidean statistics.
#' @param data Dataframe contain significant sRNA of a particular class
#' or type (ie. 24-nt or 21/22-nt sRNA)
#' @param colours The colours used to produce the heatmap image.
#' Defaults to heat colors from `grDevices` (heat.colors(100)).
#' @param dendogram logical indicating whether to include the dendrogram, and
#' retain clustering.
#' @details The function create a heatmap based on the hierarchical clustering
#' of FPKM values using euclidean statistics.
#' @examples
#'
#' ## DESeq2 example: mobile 24-nt & 21/22-nt sRNA
#'
#'  # plot heatmap of likely mobile 24-nt sRNA
#'  data("sRNA_24_mobile_DESeq2")
#'  p1 <-  heatmap_plot(sRNA_24_mobile_DESeq2)
#'
#'  # plot heatmap of likely mobile 24-nt sRNA
#'  data("sRNA_2122_mobile_DESeq2")
#'  p2 <-  heatmap_plot(sRNA_2122_mobile_DESeq2)
#'
#' ## edgeR example: mobile 24-nt & 21/22-nt sRNA
#'
#'  # plot heatmap of likely mobile 24-nt sRNA
#'  data("sRNA_24_mobile_edgeR")
#'  p3 <-  heatmap_plot(sRNA_24_mobile_edgeR)
#'
#'  # plot heatmap of likely mobile 24-nt sRNA
#'  data("sRNA_2122_mobile_edgeR")
#'  p4 <-  heatmap_plot(sRNA_2122_mobile_edgeR)
#'
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom dplyr "select"
#' @importFrom tidyselect "starts_with"
#' @importFrom stats "dist"
#' @importFrom stats "hclust"
#' @importFrom stats "as.dendrogram"
#' @importFrom stats "reorder"
#' @importFrom gplots "heatmap.2"
#' @importFrom stats "na.omit"
#' @importFrom grDevices "heat.colors"

heatmap_plot <- function(data, colours = NULL, dendogram = TRUE){
  if (base::missing(data) || !base::inherits(data, c("matrix",
                                                     "data.frame",
                                                     "DataFrame"))) {
    stop("data must be an object of class matrix, data.frame,
         DataFrame. See ?heatmap_plot for more information.")
  }
  data <- data %>% dplyr::select(!FPKM_mean)%>%
    dplyr::select(tidyselect::starts_with("FPKM"))
  matrix<- as.matrix(data)
  data[data == 0] <- 0.0001  # log trans
  data <- base::log(data,2)
  data <- stats::na.omit(data)   # remove nas.
  v <- seq(1,nrow(data), by = 1)   # subset loci clusters
  distance <- stats::dist(data[v,], method = "euclidean") # distance matrix
  cluster <- stats::hclust(distance,method="ward.D")
  # define dendrogram
  dendrogram <- stats::as.dendrogram(cluster)
  # row means
  rowv <- base::rowMeans(data, na.rm = T)
  drow <- stats::reorder(dendrogram, rowv)
  # This reorders the dendrogram as much as possible based on row/column mean.
  #This enables you to reproduce the default order with this additional step.
  #iE, allows you to order the rows of the heatmap produced
  reorderfun = function(d,w) { d }
  if (is.null(colours)) {
    plot.colours <- grDevices::heat.colors(100)
  }
  if (dendogram == TRUE) {
    p1 <- gplots::heatmap.2(as.matrix(data[v,]),
                            Rowv= dendrogram,
                            Colv=T, # t
                            dendrogram="both",
                            scale="none",
                            density.info="none",
                            trace="none",
                            reorderfun=reorderfun,
                            col= plot.colours,
                            margins =c(8,5),
                            cexCol = 1,
                            lhei = c(3,8))
  } else
    if (dendogram == FALSE) {
      p1 <- gplots::heatmap.2(as.matrix(data[v,]),
                              Rowv= TRUE, # is TRUE, which implies dendrogram is computed and reordered based on row means.
                              Colv=TRUE, #columns should be treated identically to the rows.
                              dendrogram="none",
                              scale="none",
                              density.info="none",
                              trace="none",
                              reorderfun=reorderfun,
                              col= plot.colours,
                              margins =c(8,5),
                              cexCol = 1,
                              lhei = c(3,8))
    }


  return(p1)
}








