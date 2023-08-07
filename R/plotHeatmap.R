#' Heatmap of log-transformed RPM data
#'
#' @description Undertakes RPM/FPKM normalisation using a pseudocount and 
#' transforms the data using log-scale, to enable visualization of the 
#' differences and patterns in expression across samples using a heatmap. 
#'
#'
#' @param data Dataframe; Must follow the structure created by 
#' `mobileRNA::RNAimport()`. 
#'
#' @param colours vector of colors. Default is `viridis::viridis(100)`
#'
#' @param pseudocount numeric; pseudocount, default is `1e-6`
#' 
#' @param cluster logical; include hierarchical clustering when default 
#' `cluster= TRUE`
#' 
#' @param scale character; indicating whether the values should be centered & 
#' scaled in either the row direction or the column direction, or none. 
#' Respective options are "row", "column" & "none". Default, `scale="none"`.
#' @param clustering_method Character; clustering method used. Accepts the same 
#' values as hclust. Default `clustering_method= "complete"`
#' @param row.names logical; indicated whether to include cluster names as 
#' rownames. Default `row.names=TRUE`
#'
#' @details Undertakes FPKM/RPM normalisation using a pseudocount and then 
#' transforms the normalised-RPM data using log-scale. The 
#' `mobileRNA::RNAimport()` function will  include the RPM data 
#' columns for sRNAseq data processed by ShortStack, while with mRNAseq data 
#' the function calculates FPKM.
#' The data is then plotted as a heatmap, utilising the `pheatmaps` package.  
#' 
#' This function expects to receive a dataframe containing FPKM/RPM data from 
#' sRNA-seq studies. This function employs the use of a pseudocount during 
#' normalisation as the function is expected to be used when identifying mobile 
#' sRNAs in a chimeric system. In such system, it is expected that control 
#' replicates will contain zero values for the candidate mobile sRNA clusters. 
#' 
#' 
#'@return Produces a list of objects, including the plot.
#'
#' @examples
#'
#' data("sRNA_data_mobile")
#'
#' # plot heatmap of potential mobile sRNAs
#'  p1 <-  plotHeatmap(sRNA_data_mobile)
#'
#'
#'
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
#' @importFrom viridis "viridis"
plotHeatmap <- function (data, pseudocount = 1e-6, 
                         colours = viridis::viridis(100), 
                         cluster = TRUE, scale = "none", 
                         clustering_method = "complete", 
                         row.names = TRUE) 
{
  if (base::missing(data) || !base::inherits(data, c("matrix", 
                                                     "data.frame", "DataFrame"))) {
    stop("data must be an object of class matrix, data.frame,\n DataFrame. See ?plotHeatmap for more information.")
  }
  if(any(grepl(paste0("^", "FPKM_"), colnames(data)))){
    select_data <- data %>% dplyr::select(tidyselect::starts_with("FPKM"))
    names(select_data) <- sub('^FPKM', '', names(select_data)) 
  } else 
    if(any(grepl(paste0("^", "RPM_"), colnames(data)))){
      select_data <- data %>% dplyr::select(tidyselect::starts_with("RPM_"))
      names(select_data) <- sub('^RPM_', '', names(select_data))
    } else {
      stop("data must contain columns containing either FPKM or RPM data columns.")
    }
  rownames(select_data) <- data$clusterID
  # RPM normalization with pseudocount addition
  total_reads_per_sample <- colSums(select_data)
  rpm_matrix <- (select_data / (total_reads_per_sample + pseudocount)) * 1e6
  # log transform. 
  log_rpm_matrix <- log2(rpm_matrix + 1)
  # add cluster names
  data_rownames <- data$Cluster
  # add row names 
  rownames(log_rpm_matrix) <- data_rownames
  if(cluster == FALSE){
    p1 <- pheatmap::pheatmap(log_rpm_matrix,
                             scale = scale,               
                             cluster_rows = FALSE, 
                             cluster_cols = FALSE, 
                             show_row_dendrogram = FALSE, 
                             show_col_dendrogram = FALSE,  
                             color = colours,
                             show_rownames = row.names,
                             fontsize_row = 10,            
                             fontsize_col = 10)
  } else {
    p1 <- pheatmap::pheatmap(log_rpm_matrix,
                             scale = scale,               
                             clustering_method = clustering_method,  
                             color = colours,
                             show_rownames = row.names,
                             fontsize_row = 10,            
                             fontsize_col = 10, 
                             cluster_cols = FALSE)
  }
  out <- list(plot = p1, data = log_rpm_matrix)
  return(out)
}