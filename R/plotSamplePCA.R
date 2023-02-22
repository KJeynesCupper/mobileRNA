#' PCA plot to illustrate sample distance.
#'
#' @description A function to draw a PCA plot. The function undertakes rlog
#'  transformation of the data in an unbiased manner (\code{blind=TRUE}).
#'
#' @param data class data.frame; data set containing the raw data produced as
#' output from `RNAconsensus()` and/or `RNAsubset()`.
#'
#' @param group character vector containing experimental conditions for each
#' replicate. Ensure this is in the same order as the replicates are found in
#' the `data` augment (from left to right).
#' @return A PCA plot to show sample distance.
#'
#' @details This function utilises the DESeq2 package to organise and plot the
#' data. It organises the data into a DESeqDataSet which undergoes log
#' transformation where the results are used to undertake the PCA analysis. The
#' results are plotted against the priniciple components 1 and 2.
#' @examples
#' data(sRNA_24)
#'
#' groups <- c("Tomato/Eggplant", "Tomato/Eggplant", "Tomato/Eggplant",
#'             "Tomato/Tomato",   "Tomato/Tomato", "Tomato/Tomato")
#' p <-  plotSamplePCA(data = sRNA_24,group = groups )
#'

#'
#'
#'
#' @export
#' @importFrom dplyr "select"
#' @importFrom DESeq2 "DESeqDataSetFromMatrix"
#' @importFrom SimDesign "quiet"
#' @importFrom stats "relevel"
#' @importFrom DESeq2 "estimateSizeFactors"
#' @importFrom DESeq2 "rlog"
#' @importFrom DESeq2 "plotPCA"
#' @importFrom ggplot2 "position_nudge"
#' @importFrom ggrepel "geom_label_repel"
#' @importFrom ggplot2 "aes"
plotSamplePCA <- function(data, group){
  message("Checking and organising data")
  if (base::missing(data) || !base::inherits(data, c("data.frame"))) {
    stop("data must be an object of class data.frame containing raw count data")
  }
  if (base::missing(group) || !base::inherits(group, c("character"))) {
    stop("group must be an object of class character vector containing the
         experimental condition (Treatment vs. Control)")
  }
  data <- data %>% dplyr::select(tidyselect::starts_with("Count"))
  # use DESeq to organise the data.
  column.data <- data.frame(conditions=as.factor(group))
  base::rownames(column.data) <- base::colnames(data)
  count.data.set <- SimDesign::quiet(DESeq2::DESeqDataSetFromMatrix(countData=data,
                                                   colData=column.data,
                                                   design= ~conditions))
  count.data.set$conditions <- stats::relevel(count.data.set$conditions,
                                              group[1])
  dds <- DESeq2::estimateSizeFactors(count.data.set)
  message("Transforming the count data to the log2 scale")
  # log transform the data.
  rld1 <- DESeq2::rlog(dds, blind = TRUE)
  # use the DEseq plot pca function, store in an object.
  pca <- DESeq2::plotPCA(rld1, returnData = TRUE, intgroup = "conditions")
  ## change position
  nudge <- ggplot2::position_nudge(y = 1)
  df_sample_names <- base::names(data$sizeFactor)
  sample_names <- colnames(data)
  sample_names <- sub("Count_", "", sample_names)
  pca["ID"] <- sample_names # create new column with sample names
  message("Organising principal component analysis")
  X <-DESeq2::plotPCA(rld1, intgroup = "conditions")+
  ggrepel::geom_label_repel(data = pca, ggplot2::aes(label = ID),
                            position = nudge, show.legend = FALSE)
return(X)
}

