#' Extract the distribution of RNA classes within each sample
#'
#' @description ` RNAdistribution`  extracts information from the RNA summary
#' data and calculates the total number of each RNA class identified within a
#' sample, for all samples. The results are plotted as either bar chart, with
#' an option to facet all plots into a single plot, or as a line graph.
#' @param data a dataframe, on which one of the following functions has already
#' been called: [RNAimport()],[RNAconsensus()].
#'
#' @param samples character vector. Store names of samples to analyse and plot.
#' Functionality is dependent on the use of argument `facet` for a bar chart
#' style plot and `together` for a line graph style plot.
#'
#' @param facet Logical; forms a matrix of panels defined by row and column faceting
#' variables. It plots the results for each sample as a bar-chart and contains
#' it within a single plot. The number of rows in the facet can be changed using
#' the argument ` facet.arrange` . Default ` facet = TRUE` , plot each
#' sample separately when ` facet = FALSE` .
#'
#' @param facet.arrange numeric sent to the ` ncol`  argument in
#' [facet_wrap()] to define the number of columns.
#'
#' @param colour bar plot fill colour. Default colour is "darkblue".
#'
#' @param style plotting option to choose the style of a line graph or bar chart
#' to represent your data. Where `style="line"` a line graph will be produced,
#' while `style=bar` produces a bar graph.
#'
#' @param together Logical; forms a single line graph with multiple lines each
#' to represent the sample replicates. Default `together=TRUE`.
#'
#' @return Returns a list containing the data frame and the plot(s). To access
#' one of the elements simply use the "$" symbol, and the elements "data" and
#' "plot" will appear. The `samples` argument allows uses to plot specific
#' samples in a single plot (facet bar plot or line graph). This can encourage
#' closer comparision between sample replicates.
#'
#'
#' @examples
#' data('sRNA_data')
#'
#' p1 <- RNAdistribution(data = sRNA_data, style = "line")
#'
#' p1.2 <- RNAdistribution(data = sRNA_data, style = "line",
#'                         samples = c("TomEgg_1", "TomEgg_2", "TomEgg_3"))
#' p2 <- RNAdistribution(data = sRNA_data, style = "line", together =FALSE )
#'
#' p3 <- RNAdistribution(data = sRNA_data, style = "bar")
#'
#' p3.2 <- RNAdistribution(data = sRNA_data, style = "bar",
#'                        samples = c("TomEgg_1", "TomEgg_2", "TomEgg_3"))
#'
#' p4 <- RNAdistribution(data = sRNA_data, style = "bar", facet = FALSE)
#'
#' p4 <- RNAdistribution(data = sRNA_data, style = "bar",
#'                       facet = FALSE, facet.arrange = 2 )
#'
#' @export
#' @importFrom BiocGenerics "grep"
#' @importFrom dplyr "%>%"
#' @importFrom dplyr "select"
#' @importFrom tidyselect "all_of"
#' @importFrom data.table "setDT"
#' @importFrom ggplot2 "ggplot"
#' @importFrom ggplot2 "geom_bar"
#' @importFrom ggplot2 "labs"
#' @importFrom ggplot2 "aes"
#' @importFrom ggplot2 "aes_string"
#' @importFrom tidyr "gather"
#' @importFrom data.table "melt"
#' @importFrom ggplot2 "facet_wrap"
#' @importFrom ggplot2 "theme_classic"
#' @importFrom ggplot2 "geom_point"
#' @importFrom ggplot2 "geom_line"
#' @importFrom ggplot2 "xlab"
#' @importFrom ggplot2 "ylab"
#' @importFrom ggplot2 "labs"
#'
RNAdistribution <- function(data,
                            samples =NULL,
                            style = c("bar", "line"),
                            facet = TRUE,
                            facet.arrange = 3,
                            colour = NULL, together = TRUE){

  if (base::missing(data)|| !base::inherits(data, c("data.frame"))) {
    stop("data must be a data frame, see ?help for more details")
  }
  if (base::missing(style) || !style %in% c("bar", "line")) {
    stop(paste("Please specify type of plot to produce", "(\"line\", or
               \"bar\")"))
  }

  # subset dataframe to contain specified columns
  data.cols <- data %>%
    dplyr::select(tidyselect::starts_with("DicerCall_"))
  # make list to store results
  counts.df <- apply(data.cols,MARGIN = 2,table)
  # make as data frame
  counts.df <- data.frame(counts.df)
  # remove the dicer part from names, so is only the sample names
  colnames(counts.df)<-gsub("DicerCall_","",colnames(counts.df))
  # print results
  counts.df <- data.table::setDT(counts.df, keep.rownames = "Class")[]
  # store plots in a list - can access individually
  if (is.null(colour)) {
    colour <- "darkblue"
  }

  style <- base::match.arg(style)
  if (style == "bar") {
  plist = sapply(names(counts.df)[-grep("Class", names(counts.df))],
                 function(col) {
    ggplot2::ggplot(counts.df, ggplot2::aes_string(x = "Class", y = col)) +
      ggplot2::geom_bar(stat = "identity", fill = colour)+ theme_classic()+
      ggplot2::labs(title = col, x  = "RNA Class", y = "Count") },
    simplify=FALSE)
  sn <- names(plist)
  # plot results
  if (facet == TRUE) {
    message("Printing plots as facet for samples: ", paste(sn, collapse=", "))
    # REMOVE COLUMSN IS SPECIFIC SAMPLES SPECIFIC
    if (is.null(samples)) {
      p <- print(ggplot2::ggplot(tidyr::gather(counts.df, key, Count, -Class),
                                 ggplot2::aes(Class, Count)) +
                   ggplot2::geom_bar(stat = "identity", fill = colour) +
                   ggplot2::theme_classic()+
                   ggplot2::facet_wrap(~ key, scales="free_y",
                                       ncol=facet.arrange))
    } else
      if(!is.null(samples)){
        counts.df <- counts.df %>% select(!all_of(samples))
        p <- print(ggplot2::ggplot(tidyr::gather(counts.df, key, Count, -Class),
                                   ggplot2::aes(Class, Count)) +
                     ggplot2::geom_bar(stat = "identity", fill = colour) +
                     ggplot2::theme_classic()+
                     ggplot2::facet_wrap(~ key, scales="free_y",
                                         ncol=facet.arrange))
      }

  } else
    if (facet == FALSE){ # plot individually
      message("Printing plots for samples: ", paste(sn, collapse=", "),
              domain = NULL)
     for (i in 1:length(plist)) {
       print(plist[i])
     }
    }
  }
  if (style == "line") {
    if (together == TRUE){
      if (is.null(samples)) {
        counts.df <- data.table::melt(counts.df, id.vars="Class")
        p <- print(ggplot2::ggplot(counts.df, ggplot2::aes(Class,value, col=variable, group=1)) +
                     ggplot2::geom_point() +
                     ggplot2::geom_line() +
                     ggplot2::theme_classic()+
                     ggplot2::xlab("RNA Class") +
                     ggplot2::ylab("Counts") +
                     ggplot2::labs(color='Samples'))
      } else
        if (!is.null(samples)) {
          counts.df <- counts.df %>% select(!all_of(samples))
          counts.df <- data.table::melt(counts.df, id.vars="Class")
          p <- print(ggplot2::ggplot(counts.df, ggplot2::aes(Class,value, col=variable, group=1)) +
                       ggplot2::geom_point() +
                       ggplot2::geom_line() +
                       ggplot2::theme_classic()+
                       ggplot2::xlab("RNA Class") +
                       ggplot2::ylab("Counts") +
                       ggplot2::labs(color='Samples'))
        }
    } else
    if (together == FALSE){
      plist2 = sapply(names(counts.df)[-grep("Class", names(counts.df))],
                     function(col) {
                       ggplot2::ggplot(counts.df,
                                       ggplot2::aes_string(x = "Class",
                                                           y = col, group=1)) +
                         ggplot2::geom_point(colour = colour) +
                         ggplot2::theme_classic()+
                         ggplot2::geom_line(colour = colour) +
                         ggplot2::labs(title = col,
                                       x  = "RNA Class", y = "Count")},
                     simplify=FALSE)
      sn <- names(plist2)
      p <-plist2
      message("Printing line plots for samples: ", paste(sn, collapse=", "))
      for (J in 1:length(plist2)) {
        print(plist2[J])
      }
    }
  }
  out <- list(plot = p, data = counts.df)
  print(counts.df)
  return(out)
}


