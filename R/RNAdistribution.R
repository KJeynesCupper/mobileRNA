#' Plot the distribution of sRNA lengths
#'
#' @description \code{RNAdistribution} plots the distribution of dicer-derived
#' sRNA classes (20-24nt) across samples or across the sRNA consensus
#' determined by the function [mobileRNA::RNAconsensus()].
#'
#' @param data a dataframe, on which one of the following functions has already
#' been called: [mobileRNA::RNAimport()],[mobileRNA::RNAconsensus()].
#'
#' @param samples character vector. Store names of samples to analyse and plot.
#' Argument is required for plotting individual sample replicates, either
#' individually, overlapped together or in a facet. Use the sample
#' replicate names present in the data frame (`data`), select samples you wish
#' to plot. Is not required when plotting the sRNA consensus using the argument
#' \code{total=TRUE}.
#'
#' @param facet Logical; forms a matrix of panels defined by row and column faceting
#' variables. It plots the results for each sample as a bar-chart and contains
#' it within a single plot. The number of rows in the facet can be changed using
#' the argument \code{facet.arrange} . Default \code{facet = TRUE} , plot each
#' sample separately when \code{facet = FALSE} .
#'
#' @param facet.arrange numeric; value supplied to define the number of columns
#' to include in the facet. This argument is piped into the ` ncol`  argument in
#' [ggplot2::facet_wrap()] to define the number of columns. By default, this is
#' set to 3.
#'
#' @param colour bar plot fill colour. Default colour is "darkblue".
#'
#' @param style plotting option to choose the style of either a line graph or
#' bar chart to represent your data.
#' * Where \code{style="line"} a line graph will be produced
#' * Where \code{style="bar"} produces a bar graph
#' * Where \code{style="consensus"} produces the line graph for the consenus
#' sRNA in conjunction with \code{consensus=TRUE}
#'
#' @param together Logical; forms a single line graph with multiple lines each
#' to represent the sample replicates. Default \code{together=TRUE}.
#'
#'
#'@param consensus Logical; plots the distribution of sRNA classes across all
#'identified dicer-derived clusters based on the consensus. See
#'[mobileRNA::RNAconsensus()]
#'function to calculate consensus RNA class based the experimental replicates.
#'Default \code{consensus=FALSE}.
#'
#'@param relative Logical; calculates relative frequency of consensus
#'dicer-derived sRNA clusters. Only applicable when  only in conjunction with
#' argument `consensus`, when set as \code{consensus=TRUE}. As default,
#' \code{relative=FALSE}.
#'
#'@details
#'
#'The function can be used to plot a variety of different comparisons and plots.
#'It can be used to plot the distribution of sRNA classes within each sample
#'replicate, which can be represented as a bar chart \code{style="bar"} or a
#'line graph \code{style="line"}. These plots can be represented individually or
#'in a single plot facet \code{facet="TRUE"} by default.
#'
#'
#'To plot the sRNA dicer-derived clusters identified in each sample, the
#'function extracts the information from the RNA summary data and calculates the
#'total number of each RNA class identified within a sample, for all samples.
#'
#'Alternatively, the function allows you to plot the line graph for each sample
#'together, overlapped on a single graph \code{total="TRUE"}. This is not an
#'option for bar plots.
#'
#'The final option, is to plot the total consensus of dicer-derived sRNA clusters
#'across the experimental conditions, the function pulls the consensus call
#'from the column created by the [RNAconsensus()] function in the working
#'data frame.
#'
#'
#'
#'@return The function returns a list containing the results in the form of a
#'data frame and the plot(s). To access an element, simply use the "$" symbol,
#'and the elements "data" and "plot" will appear. The `samples` argument allows
#'uses to plot specific samples in a single plot (facet bar plot or line graph).
#'This can encourage closer comparison between sample replicates.
#'
#' @examples
#' data('sRNA_data')
#'
#' p1 <- RNAdistribution(data = sRNA_data, style = "line")
#'
#' p1.2 <- RNAdistribution(data = sRNA_data, style = "line",
#'                         samples = c("heterograft_1", "heterograft_2", "heterograft_3"))
#' p2 <- RNAdistribution(data = sRNA_data, style = "line", together =FALSE )
#'
#' p3 <- RNAdistribution(data = sRNA_data, style = "bar")
#'
#' p3.2 <- RNAdistribution(data = sRNA_data, style = "bar",
#'                        samples = c("heterograft_1", "heterograft_2", "heterograft_3"))
#'
#' p4 <- RNAdistribution(data = sRNA_data, style = "bar", facet = FALSE)
#'
#' p5 <- RNAdistribution(data = sRNA_data, style = "bar",
#'                       facet = TRUE, facet.arrange = 2 )
#'
#'data("sRNA_data_consensus")
#'p6 <- RNAdistribution(data = sRNA_data_consensus, style = "consensus", consensus = TRUE)
#'
#' @export
#' @importFrom BiocGenerics "grep"
#' @importFrom dplyr "%>%"
#' @importFrom dplyr "select"
#' @importFrom dplyr "count"
#' @importFrom dplyr "mutate"
#' @importFrom tidyselect "all_of"
#' @importFrom data.table "setDT"
#' @importFrom ggplot2 "ggplot"
#' @importFrom ggplot2 "aes_string"
#' @importFrom ggplot2 "geom_bar"
#' @importFrom ggplot2 "labs"
#' @importFrom ggplot2 "aes"
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
RNAdistribution <- RNAdistribution2 <- function(data, samples =NULL,
                                                style = c("bar", "line", "consensus"),
                                                facet = TRUE,
                                                facet.arrange = 3,
                                                colour = "darkblue", together = TRUE,
                                                consensus = FALSE, relative = FALSE) {

  if (base::missing(data)|| !base::inherits(data, c("data.frame"))) {
    stop("data must be a data frame, see ?help for more details")
  }
  if (base::missing(style) || !style %in% c("bar", "line", "consensus")) {
    stop(paste("Please specify type of plot to produce", "(\"line\", or
               \"bar\")"))
  }
  # PLOT CONSENSUS
  if (consensus == TRUE){
    x <- data %>% dplyr::count(sRNA_Consensus)
    if(!relative == FALSE){
      x <- x %>% dplyr::mutate(freq = n / sum(n))
      p1 <- ggplot2::ggplot(x, ggplot2::aes(x = sRNA_Consensus, y = freq, group=1)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::theme_classic()+
        ggplot2::xlab("RNA Class") +
        ggplot2::ylab("Relative frequency")
    } else
      p1 <- ggplot2::ggplot(x, ggplot2::aes(x = sRNA_Consensus, y = n, group=1)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::theme_classic()+
        ggplot2::xlab("RNA Class") +
        ggplot2::ylab("Counts")

    out <- list(plot = p1, data = x)
    return(out)
    print(x)
    print(p1)
  } else {

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
    # Remove row with unclassified sRNA (class = N)
    counts.df<-counts.df[!(counts.df$Class=="N"),]

    style <- base::match.arg(style)
    if (style == "bar") {
      plist = sapply(names(counts.df)[-grep("Class", names(counts.df))],
                     function(col) {
                       ggplot2::ggplot(counts.df, ggplot2::aes_string(x = "Class", y = col)) +
                         ggplot2::geom_bar(stat = "identity", fill = colour)+
                         ggplot2::theme_classic()+
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
        print(counts.df)
        out <- list(plot = p, data = counts.df)
        return(out)

      } else
        if (facet == FALSE){ # plot individually
          message("Printing plots for samples: ", paste(sn, collapse=", "),
                  domain = NULL)
          save <- list()
          for (i in 1:length(plist)) {
            print(plist[i])
            save <- list(save, plist[i])
          }
          print(counts.df)
          out <- list(plot = save, data = counts.df)
          return(out)
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
            out <- list(plot = p, data = counts.df)
            return(out)
          }
      } else
        if (together == FALSE){
          if(facet == TRUE){
            message("Printing plots as facet for samples")
            if (is.null(samples)) {
              p <- print(ggplot2::ggplot(tidyr::gather(counts.df, key, Count, -Class),
                                         ggplot2::aes(Class, Count, group=1)) +
                           ggplot2::geom_point() +
                           ggplot2::geom_line() +
                           ggplot2::theme_classic()+
                           ggplot2::facet_wrap(~ key, scales="free_y",
                                               ncol=facet.arrange))
            } else
              if(!is.null(samples)){
                counts.df <- counts.df %>% select(!all_of(samples))
                p <- print(ggplot2::ggplot(tidyr::gather(counts.df, key, Count, -Class),
                                           ggplot2::aes(Class, Count, group=1)) +
                             ggplot2::geom_point() +
                             ggplot2::geom_line() +
                             ggplot2::theme_classic()+
                             ggplot2::facet_wrap(~ key, scales="free_y",
                                                 ncol=facet.arrange))
              }
            #print(counts.df)
            out <- list(plot = p, data = counts.df)
            return(out)
          }
          else
            if (facet == FALSE){
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
              out <- list(plot = p, data = counts.df)
              return(out)

              message("Printing line plots for samples: ", paste(sn, collapse=", "))
              for (J in 1:length(plist2)) {
                print(plist2[J])
              }
            }
        }
    }
  }
}
