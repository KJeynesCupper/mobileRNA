#' Volcano plots (logFC, FDR)
#'
#' @description Volcano plots shows the relationship between fold change
#' and adjusted p-value (i.e. FDR). Each point represents an individual genes or 
#' cluster. Pointes are coloured according to their significance level 
#' thresholds (0.05, 0.01, and 0.001 or unchanged). 
#'
#' @param data Dataframe; Must follow the structure created by 
#' `mobileRNA::RNAimport()` and contain columns contributed by statistical 
#' analysis using `mobileRNA::RNAanalysis()`
#'
#' @param labels Add labels molecules with the most significant changes
#'
#' @param top.molecules numeric; number of molecules to label
#' 
#' @param colour.scheme colour scale; recommend picking a 
#' "Viridis colour scales from viridisLite". Default set to 
#' ggplot2::scale_colour_viridis_d
#' 
#' 
#' @export
#' @importFrom dplyr "mutate"
#' @importFrom dplyr "case_when"
#' @importFrom ggplot2 "ggplot"
#' @importFrom ggplot2 "aes"
#' @importFrom ggplot2 "geom_point"
#' @importFrom ggplot2 "xlab"
#' @importFrom ggplot2 "ylab"
#' @importFrom ggplot2 "scale_color_viridis_d"
#' @importFrom ggplot2 "guides"
#' @importFrom ggplot2 "guide_legend"
#' @importFrom dplyr "filter"
#' @importFrom dplyr "arrange"
#' @importFrom dplyr "desc"
#' @importFrom dplyr "count"
#' @importFrom utils "head"
#' @importFrom dplyr "bind_rows"
#' @importFrom ggplot2 "scale_color_viridis_d"
plotVolcano <- function(data, labels = FALSE, top.molecules = 10, 
                        colour.scheme = NULL){
  if (base::missing(data) || !base::inherits(data, c("data.frame"))) {
    stop("data must be an object of class data.frame containing raw count data")
  }
  
  # add expression column to state up or down, ass significance level col. 
  data_val <- data %>% 
    dplyr::mutate(
      expression = dplyr:: case_when(log2FoldChange >= log(2) & padjusted <=0.05 
                                     ~ "Up-regulated",
                             log2FoldChange <= -log(2) & padjusted <= 0.05
                             ~ "Down-regulated",
                             TRUE ~ "Unchanged")) %>% 
    dplyr::mutate(
          significance = dplyr::case_when(
            abs(log2FoldChange) >= log(2) & padjusted <= 0.05 & padjusted > 0.01 
            ~ "FDR 0.05", 
            abs(log2FoldChange) >= log(2) & padjusted <= 0.01 & padjusted >0.001
            ~ "FDR 0.01",
            abs(log2FoldChange) >= log(2) & padjusted <= 0.001 ~ "FDR 0.001", 
            TRUE ~ "Unchanged"))
      
  if(is.null(colour.scheme)){
        colour.scheme <- ggplot2::scale_color_viridis_d()
      }
 plot_out <- ggplot2::ggplot(data_val, ggplot2::aes(log2FoldChange, 
                                                    -log(padjusted,10))) +
  ggplot2::geom_point(ggplot2::aes(color = significance), size = 2/5)+
  ggplot2::xlab(expression("log"[2]*"FC"))+
  ggplot2::ylab(expression("-log"[10]*"FDR"))+
  colour.scheme +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=1.5)))

if(labels){
  top_molecules <- dplyr::bind_rows(
    data_val %>% 
      dplyr::filter(expression == 'Upregulated') %>% 
      dplyr::arrange(padjusted, dplyr::desc(abs(log2FoldChange))) %>% 
      utils::head(top.molecules),
    data_val %>% 
      dplyr::filter(expression == 'Downregulated') %>% 
      dplyr::arrange(padjusted, dplyr::desc(abs(log2FoldChange))) %>% 
      utils::head(top.molecules)
  )
  
  plot_out <-  plot_out +
    ggrepel::geom_label_repel(data = top_molecules,
                     mapping = ggplot2::aes(log2FoldChange, -log(padjusted,10), 
                                            label = Locus),size = 2)
}
      
data_summary<-  data_val %>% 
  dplyr::count(expression, significance) 

out <- list(plot = plot_out, data = data_summary)
return(out)
}