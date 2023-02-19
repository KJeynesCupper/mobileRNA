#' Plot Consensus small RNA frequency
#'
#' @description Plot a bar chart of distribution of small RNA classes
#' (20-24nt sRNA) identified in the samples.
#'
#' @param data A data-frame (See [sample_table()] and
#' [define_consensus()] to produce an appropriate data-frame).
#' Plots data based on column labeled `"sRNA_Consensus"`.
#' @param relative If `relative=TRUE`, will plot the relative frequency.
#' As default, the function will plot the absolute frequency
#' (`relative == FALSE`) rather than the relative frequency.
#' @return A bar chart of the distribution of consensus sRNA in the data.
#' @examples
#'
#' data("sRNA_data_summary")
#' p1 <- frequency_plot(data = sRNA_data_summary, relative=TRUE )
#'
#'
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom dplyr "count"
#' @importFrom dplyr "mutate"
#' @importFrom ggplot2 "ggplot"
#' @importFrom ggplot2 "geom_bar"
#' @importFrom ggplot2 "labs"
#' @importFrom ggplot2 "theme_classic"
frequency_plot <- function(data, relative = FALSE){
  if (base::missing(data) || !base::inherits(data, c("matrix", "data.frame",
                                                     "DataFrame"))) {
    stop("data must be an object of class matrix, data.frame, DataFrame.
         See ?frequency_plot for more information")
  }
  x <- data %>% dplyr::count(sRNA_Consensus)
  if(!relative == FALSE){
    x <- x %>% dplyr::mutate(freq = n / sum(n))
  }
  p1 <- ggplot2::ggplot(data= x, ggplot2::aes(x = sRNA_Consensus,
                                              y = n)) +
    ggplot2::geom_bar(stat = "identity", fill = "black")+
    ggplot2::labs(x  = "siRNA class", y = "Occurances",
                  title = "The frequency of siRNA classes")+
    ggplot2::theme_classic()
  return(p1)
}







