#' Subset data to select a specific small RNA population
#'
#' @description Creates a data-frame containing the desired sRNA class(s) based
#' on the consensus sRNA determination.
#'
#' @details
#' See [mobileRNA::RNAdicercall()] for information on defining the sRNA type
#' for each cluster. The function allows the choice to filtered the data by
#' statistical significance based on differential expression analysis, see
#' [mobileRNA::RNAanalysis()]. Set \code{sig=TRUE} to filtered by significance
#' (p-adjusted). It is important to consider the point in your analysis you
#' subset the data or/and undertake differential analysis to achieve statistical
#' values. Subsetting the dataset into groups based on the sRNA class will
#' create a smaller set of data for each to draw statistical differences.
#' Depending on the size of your data, and analysis aims this should be taken
#' into consideration.
#'
#'
#'
#'
#' @param data A numerical data-frame containing the sample data, with a
#' defined consensus sRNA class/type for each sRNA dicer-derived cluster
#' (see [mobileRNA::RNAdicercall()].
#'
#' @param type numeric; small RNA class(es) to select.
#' 
#' @param sig Parameter to filter and select significant sRNA. If
#'  \code{sig=TRUE}, data will be filtered based on p-adjusted < 0.05
#'  significance threshold.
#'
#' @return A dataframe containing sRNA clusters with a sRNA consensus matching
#' the size instructed to the type argument.
#' @examples
#' data("sRNA_data_consensus")
#'
#' # Subset data for  24-nt sRNAs
#' sRNA_24 <- RNAsubset(sRNA_data_consensus, type = 24)
#'
#'
#'# Subset data for 24 21/22-nt sRNAs
#'sRNA_2122 <- RNAsubset(sRNA_data_consensus, type = c(21, 22))
#'
#' # You can subset by any combination of classes. For example, a dataset
#' # of 23-nt & 24-nt sRNAs or just 20-nt sRNAs.
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr "filter"

RNAsubset <- function(data, type,  sig=FALSE){
    x <- data %>% dplyr::filter(DicerConsensus %in% type)
    if(sig){
      x<- x %>%
        dplyr::filter(padj < 0.05)
    }
    return(x)
  }

