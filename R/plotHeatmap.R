#' Heatmap using hierarchical clustering
#'
#' @description Plots a heatmap with hierarchical clustering via an rlog
#' transformation of RPM data and euclidean statistics.
#'
#'
#' @param data Data frame containing FPKM values for each samples on
#' each sRNA dicer-derived cluster of interest.
#'
#'
#' @param colours  Colours to display and represent the heatmap.
#' Defaults to [grDevices::heat.colors()] (heat.colors(100)).
#'
#'
#' @param dendogram Logical; indicating whether to include the dendrogram and
#' clustering, and retain clustering. Default, \code{dendogram = TRUE} to
#' include.
#'
#'
#' @details The function create a heatmap based on the hierarchical clustering
#' of FPKM values using euclidean statistics.
#'
#'@return
#'Produces a list of objects, including the plot.
#'
#' @examples
#'
#' ## DESeq2 example: mobile 24-nt & 21/22-nt sRNA
#'
#'  # plot heatmap of likely mobile 24-nt sRNA
#'  data("sRNA_24_mobile_DESeq2")
#'  p1 <-  plotHeatmap(sRNA_24_mobile_DESeq2)
#'
#'  # plot heatmap of likely mobile 24-nt sRNA
#'  data("sRNA_2122_mobile_DESeq2")
#'  p2 <-  plotHeatmap(sRNA_2122_mobile_DESeq2)
#'
#'
#'
#' ## edgeR example: mobile 24-nt & 21/22-nt sRNA
#'
#'  # plot heatmap of likely mobile 24-nt sRNA
#'  data("sRNA_24_mobile_edgeR")
#'  p3 <-  plotHeatmap(sRNA_24_mobile_edgeR)
#'
#'  # plot heatmap of likely mobile 24-nt sRNA
#'  data("sRNA_2122_mobile_edgeR")
#'  p4 <-  plotHeatmap(sRNA_2122_mobile_edgeR)
#'
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom dplyr "select"
#' @importFrom tidyselect "starts_with"
#' @importFrom stats "dist"
#' @importFrom stats "hclust"
#' @importFrom stats "as.dendrogram"
#' @importFrom stats "reorder"
#' @importFrom pheatmap "pheatmap"
#' @importFrom stats "na.omit"
#' @importFrom grDevices "heat.colors"

plotHeatmap <-function (data, colours = NULL, dendogram = TRUE){
  if (base::missing(data) || !base::inherits(data, c("matrix",
                                                     "data.frame",
                                                     "DataFrame"))) {
    stop("data must be an object of class matrix, data.frame,\n DataFrame.
         See ?plotHeatmap for more information.")
  }
  select_data <- data %>% dplyr::select(tidyselect::starts_with("RPM_"))
  rownames(select_data) <- data$clusterID
  # remove RPM_
  for ( col in 1:ncol(select_data)){
    colnames(select_data)[col] <-  sub("RPM_", "", colnames(select_data)[col])
  }
  matrix <- as.matrix(select_data)
  select_data[select_data == 0] <- 1e-04
  select_data <- base::log(select_data, 2)
  select_data <- stats::na.omit(select_data)
  v <- seq(1, nrow(select_data), by = 1)
  distance <- stats::dist(select_data[v, ], method = "euclidean")
  cluster <- stats::hclust(distance, method = "ward.D")
  dendrogram <- stats::as.dendrogram(cluster)
  rowv <- base::rowMeans(select_data, na.rm = T)
  drow <- stats::reorder(dendrogram, rowv)
  reorderfun = function(d, w) {
    d
  }
  if (is.null(colours)) {
    plot.colours <- grDevices::heat.colors(100)
  }
  else {
    plot.colours <- colours
  }

  if (dendogram == TRUE) {

   p1 <-  pheatmap::pheatmap(as.matrix(select_data[v, ]),
                       border_color = FALSE,cellheight = 10,
                       color = plot.colours)

  }
  else if (dendogram == FALSE) {
    p2 <-  pheatmap::pheatmap(as.matrix(select_data[v, ]),
                              border_color = FALSE,
                              treeheight_row = 0,
                              treeheight_col = 0 ,
                              cluster_rows=FALSE, cluster_cols=FALSE,
                              color = plot.colours)
  }

  out <- list(plot = p1, data = as.matrix(select_data[v, ]))
  return(out)
}

