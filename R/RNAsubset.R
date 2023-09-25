#' Subset data to select a specific small RNA population
#'
#' @description Creates a data-frame containing the desired sRNA class(s) based
#' on the consensus sRNA determination.
#'
#' @details
#' See [mobileRNA::RNAdicercall()] for information on defining the sRNA class
#' for each cluster. The function allows the choice to filtered the data by
#' statistical statisticalnificance based on differential expression analysis, see
#' [mobileRNA::RNAdifferentialAnalysis()]. Set \code{statistical=TRUE} to 
#' filtered by statisticalnificance (p-adjusted). It is important to consider 
#' the point in your analysis you subset the data or/and undertake differential 
#' analysis to achieve statistical values. 
#' Subsetting the dataset into groups based on the sRNA class will
#' create a smaller set of data for each to draw statistical differences.
#' Depending on the size of your data, and analysis aims this should be taken
#' into consideration.
#'
#'
#'
#'
#' @param data A numerical data-frame containing the sample data, with a
#' defined consensus sRNA class/class for each sRNA dicer-derived cluster
#' (see [mobileRNA::RNAdicercall()].
#'
#' @param class numeric; small RNA class(es) to select.
#' 
#' @param statistical Parameter to filter and select statisticalnificant sRNA. If
#'  \code{statistical=TRUE}, data will be filtered based on p-adjusted < 0.05
#'  statisticalnificance threshold.
#'
#' @return A dataframe containing sRNA clusters with a sRNA consensus matching
#' the size instructed to the class argument.
#' @examples
#' data("sRNA_data_consensus")
#'
#' # Subset data for  24-nt sRNAs
#' sRNA_24 <- RNAsubset(sRNA_data_consensus, class = 24)
#'
#'
#'# Subset data for 24 21/22-nt sRNAs
#'sRNA_2122 <- RNAsubset(sRNA_data_consensus, class = c(21, 22))
#'
#' # You can subset by any combination of classes. For example, a dataset
#' # of 23-nt & 24-nt sRNAs or just 20-nt sRNAs.
#'
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom dplyr "filter"

RNAsubset <- function(data, class,  statistical=FALSE){
  if (!base::inherits(data, c("matrix", "data.frame", "DataFrame"))) {
    stop("data must be an object of class matrix, data.frame, DataFrame")
  }
  if (base::missing(class) || !base::inherits(class, "numeric")) {
    stop("Please specify a numeric vector stating the class of sRNAs to select")
  }
    x <- data %>% dplyr::filter(DicerConsensus %in% class)
    if(statistical){
      x<- x %>%
        dplyr::filter(padj < 0.05)
    }
    return(x)
  }

