#' Extract the distribution of RNA classes within each sample
#'
#' @description ` RNAdistribution`  extracts information from the RNA summary
#' data and calculates the total number of each RNA class identified within a
#' sample, for all samples. The results are plotted as either bar chart, with
#' an option to facet all plots into a single plot, or as a line graph.
#' @param data a dataframe, on which one of the following functions has already
#' been called: [RNAimport()],[RNAconsensus()].
#'
#' @param facet Logical; forms a matrix of panels defined by row and column faceting
#' variables. It plots the results for each sample as a bar-chart and contains
#' it within a single plot. The number of rows in the facet can be changed using
#' the augment ` facet.arrange` . Default ` facet = TRUE` , plot each
#' sample separately when ` facet = FALSE` .
#'
#' @param facet.arrange numeric sent to the ` ncol`  augment in
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
#' @examples
#' data('sRNA_data')
#' sample.distribution <- RNAdistribution(data = sRNA_data, style = "line")
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
RNAdistribution <- function(data, style = c("bar", "line"),facet = TRUE,
                            facet.arrange = 3,
                            colour = NULL, together = TRUE){
  if (base::missing(data)|| !base::inherits(data, c("data.frame"))) {
    stop("data must be a data frame, see ?help for more details")
  }
  if (base::missing(style) || !style %in% c("bar", "line")) {
    stop(paste("Please specify type of plot to produce", "(\"line\", or
               \"bar\")"))
  }
  #select columns
  class_colnames <- BiocGenerics::grep('DicerCall_', names(data), value = TRUE)
  samples <- base::sub('DicerCall_', '', class_colnames)
  # subset dataframe to contain specified columns
  data.cols <- data %>%
    dplyr::select(tidyselect::all_of(class_colnames))
  # make list to store results
  counts.df <- apply(data.cols,MARGIN = 2,table)
  # make as data frame
  counts.df <- data.frame(counts.df)
  # remove the dicer part from names, so is only the sample names
  colnames(counts.df)<-gsub("DicerCall_","",colnames(counts.df))
  # print results
  counts.df <- data.table::setDT(counts.df, keep.rownames = "Class")[]

  print(counts.df)

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
    print(ggplot2::ggplot(tidyr::gather(counts.df, key, Count, -Class),
                          ggplot2::aes(Class, Count)) +
            ggplot2::geom_bar(stat = "identity", fill = colour) +
            theme_classic()+
            ggplot2::facet_wrap(~ key, scales="free_y", ncol=facet.arrange))
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
    counts.df <- melt(counts.df, id.vars="Class")
 print(ggplot(counts.df, aes(Class,value, col=variable, group=1)) +
         geom_point() +
         geom_line() +
        theme_classic()+
        xlab("RNA Class") +
   ylab("Counts") +
   labs(color='Samples'))
    } else
    if (together == FALSE){
      plist2 = sapply(names(counts.df)[-grep("Class", names(counts.df))],
                     function(col) {
                       ggplot2::ggplot(counts.df,
                                       ggplot2::aes_string(x = "Class",
                                                           y = col, group=1)) +
                         ggplot2::geom_point(colour = colour) +
                         theme_classic()+
                         ggplot2::geom_line(colour = colour) +
                         ggplot2::labs(title = col,
                                       x  = "RNA Class", y = "Count")},
                     simplify=FALSE)
      sn <- names(plist)
      message("Printing line plots for samples: ", paste(sn, collapse=", "))
      for (J in 1:length(plist2)) {
        print(plist2[J])
      }
    }
  }
  return(counts.df)
}


