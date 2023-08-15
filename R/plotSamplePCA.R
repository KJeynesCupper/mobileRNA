#' PCA plot to illustrate sample distance.
#'
#' @description A function to draw a principal component analysis (PCA) plot.
#'  The function undertakes rlog transformation of the data in
#'  an unbiased manner (\code{blind=TRUE}).
#'
#' @param data class data.frame; data set containing the raw data produced as
#' output from [mobileRNA::RNAdicercall()] and/or [mobileRNA::RNAsubset()].
#'
#' @param group character vector; contains experimental conditions for each
#' replicate. IMPORTANT: Ensure this is in the same order as the replicates are
#' found in the data frame supplied to the `data` argument (from left to right).
#'
#' @param vst logical; Variance stabilizing transformation. By default, the
#' function uses a regularized log transformation on the data set, however, this
#' will not suit all experimental designs.
#' 
#' @param labels logical; include sample name labels on PCA. Default 
#' `labels=TRUE`
#' 
#' @param boxed logical; add box around each sample name label. Default 
#' `boxed=TRUE`
#' 
#' @param legend.title character; title for legend key. Default
#' `legend.title = "Conditions"`
#' @param size.ratio numeric; set plot ratio, broadens axis dimensions by ratio.
#' Default `size.ratio=2`, double the plot dimension. 
#'
#'@param colours Vector of HEX colour codes. Must match the number of 
#'conditions. For example, 
#'`colours = c("#E69F00", "#56B4E9", "#CC79A7", "#009E73")`
#'
#'@param point.shape Logical; set whether the point shapes should be different
#'for each condition. 
#'
#'@param ggplot.theme character; state the `ggplot2` theme (without () 
#'brackets). For example, `ggplot.theme=theme_classic`. 
#'
#' @return A PCA plot to show sample distance.
#'
#' @details This function utilises the DESeq2 package to organise and plot the
#' data. It organises the data into a DESeqDataSet which undergoes log
#' transformation where the results are used to undertake the PCA analysis. The
#' results are plotted against the principle components 1 and 2.
#'
#' In special conditions, regularized log transformation will not suit
#' the experimental design. For example, an experimental design without
#' replicates. In this instance, it is preferable to change the default setting
#' and switch to a variance stabilizing transformation method
#' (\code{`vst=TRUE`}).
#'
#'
#' @examples
#' data("sRNA_data_consensus")
#'
#' groups <- c("Heterograft", "Heterograft", "Heterograft",
#'             "Selfgraft", "Selfgraft", "Selfgraft")
#'             
#' p <-  plotSamplePCA(data = sRNA_data_consensus,group = groups )
#'
#' plot(p)
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
#' @importFrom ggplot2 "labs"
#' @importFrom ggplot2 "coord_fixed"
#' @importFrom ggrepel "geom_text_repel"
#' @importFrom ggplot2 "ggplot"
#' @importFrom ggplot2 "geom_point"
#' @importFrom ggplot2 "scale_color_manual"
plotSamplePCA <- function(data, group, vst = FALSE, labels = TRUE, boxed = TRUE,
                          legend.title = "Conditions", size.ratio = 2, 
                          colours = NULL, point.shape = TRUE, 
                          ggplot.theme = NULL){
  cat("Checking and organising data \n")
  
  if (base::missing(data) || !base::inherits(data, c("data.frame"))) {
    stop("data must be an object of class data.frame containing raw count data")
  }
  if (base::missing(group) || !base::inherits(group, c("character"))) {
    stop("group must be an object of class character vector containing the experimental condition (Treatment vs. Control)")
  }
  data <- data %>% dplyr::select(tidyselect::starts_with("Count"))
  # use DESeq to organise the data.
  column.data <- data.frame(conditions=as.factor(group))
  base::rownames(column.data) <- base::colnames(data)
  count.data.set <- SimDesign::quiet(DESeq2::DESeqDataSetFromMatrix(
    countData=data,colData=column.data,design= ~conditions))
  count.data.set$conditions <- stats::relevel(count.data.set$conditions,
                                              group[1])
  dds <- DESeq2::estimateSizeFactors(count.data.set)
  
  # log transform the data.
  
  if(vst ==TRUE){
    cat("Transforming the count data with a variance stabilizing transformation \n")
    rld1 <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
    # log transform the data.
  } else
    if(vst == FALSE) {
      cat("Transforming the count data to the log2 scale \n")
      rld1 <- DESeq2::rlog(dds, blind = TRUE) # log transform the data.
    }
  
  # use the DEseq plot pca function, store in an object.
  pca <- DESeq2::plotPCA(rld1, returnData = TRUE, intgroup = "conditions")
  ## change position
  sample_names <- sub("Count_", "", colnames(data))
  pca["ID"] <- sample_names # create new column with sample names
  percentVar <- round(100 * attr(pca, "percentVar"))
  
  cat("Organising principal component analysis \n")
  if(labels == TRUE){
    if(boxed == TRUE){
      X <- ggplot2::ggplot(pca, ggplot2::aes(PC1, PC2, color=conditions)) +
        {if(point.shape) ggplot2::geom_point(ggplot2::aes(shape = conditions, 
                                                          size=3))}+
        {if(point.shape == FALSE) ggplot2::geom_point(size=3)}+
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        {if(!is.null(colours)) ggplot2::scale_color_manual(values=colours)}+ 
        ggplot2::coord_fixed()+
        ggrepel::geom_label_repel(data = pca, ggplot2::aes(label = ID), 
                                  show.legend = FALSE, box.padding = 1)+
        ggplot2::labs(color = legend.title) + 
        ggplot2::coord_fixed(ratio = size.ratio)+
        {if(!is.null(ggplot.theme)) ggplot2::ggplot.theme() }
      
    } else
      X <- ggplot2::ggplot(pca, ggplot2::aes(PC1, PC2, color=conditions)) +
        {if(point.shape) ggplot2::geom_point(ggplot2::aes(shape = conditions, 
                                                          size=3))}+
        {if(point.shape == FALSE) ggplot2::geom_point(size=3)}+
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        {if(!is.null(colours)) ggplot2::scale_color_manual(values=colours)}+ 
        ggrepel::geom_label_repel(data = pca, ggplot2::aes(label = ID), 
                                  show.legend = FALSE, box.padding = 1)+
        ggplot2::labs(color = legend.title) + 
        suppressMessages(ggplot2::coord_fixed(ratio = size.ratio))+
        {if(!is.null(ggplot.theme)) ggplot2::ggplot.theme() }
    
    
  } else {
    X <- ggplot2::ggplot(pca, ggplot2::aes(PC1, PC2, color=conditions)) +
      {if(point.shape) ggplot2::geom_point(ggplot2::aes(shape = conditions,
                                                        size=3))}+
      {if(point.shape == FALSE) ggplot2::geom_point(size=3)}+
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      {if(!is.null(colours)) ggplot2::scale_color_manual(values=colours)}+ 
      ggplot2::labs(color = legend.title) + 
      ggplot2::coord_fixed(ratio = size.ratio)+
      {if(!is.null(ggplot.theme)) ggplot2::ggplot.theme() }
  }
  return(X)
}