#' Plot the distribution of sRNA classes based on nucleotide length
#'
#' @description \code{RNAdistribution} plots the distribution of dicer-derived
#' sRNA classes across samples or the sRNA consensus determined by the
#'  [mobileRNA::RNAdicercall()] function. This can be displayed as a line or bar 
#'  plot. 
#'
#' @param data data.frame; generated originally by [mobileRNA::RNAimport()]
#'
#' @param samples character; states a subset of samples to plot. This argument 
#' is based on the sample names within the data. 
#' 
#' @param outline character; states the outline colour for a box plot. 
#' Default is "black".
#'
#' @param facet logical; forms a matrix of panels defined by row and column
#' faceting variables. It plots the results for each sample as a bar chart
#' and contains it within a single plot. The number of rows in the facet can
#' be changed using the argument \code{facet.arrange}.
#' Default \code{facet = TRUE} , plots each sample separately when
#' \code{facet = FALSE} .
#'
#' @param facet.arrange numeric; value supplied to define the number of columns
#' to include in the facet. This argument is piped into the ` ncol`  argument in
#' [ggplot2::facet_wrap()] to define the number of columns. By default, this is
#' set to 3.
#'
#' @param colour character; fill colour, Default is "#0868AC".
#'
#'@param data.type character; either plotting "samples" or the "consensus" 
#'stored in the 'DicerConsensus' column. 
#'
#' @param style character; \code{style="line"} or \code{style="bar"}. 
#' Instructs how to plot the data. Where \code{style="line"} plots a line graph 
#' and \code{style="bar"} plots a bar graph. 
#'
#' @param overlap logical; generates a single line graph, containing all 
#' sample replicate information. Default \code{overlap=TRUE}.
#'
#'@param non.classified logical; to include distribution of sRNAs which were 
#'unclassified. Default is TRUE, ie. maintain.
#'
#'
#'@param relative logical; used in conjunction with \code{data.type="consensus"}. 
#'Instructs plotting of the relative frequency of all sRNA classes across the 
#'data defined by the consensus dicer-derived sRNA column (see 
#'[mobileRNA::RNAdicercall()] function for more information). 
#'
#'@param wrap.scales character; scales be fixed ("fixed", the default), free 
#'("free"), or free in one dimension ("free_x", "free_y")
#'
#'
#'@details
#'
#'The function can be used to plot a variety of different comparisons and plots.
#'It can be used to plot the distribution of sRNA classes within each sample
#'replicate, which can be represented as a bar chart \code{style="bar"} or a
#'line graph \code{style="line"}. These plots can be represented individually or
#'in a single facet plot when \code{facet="TRUE"}. 
#'
#'Additionally, there is the option to plot the distribution of sRNA classes
#'within individual samples or to plot the distribution of the consensus 
#' dicer-derived sRNA classes determined by the [RNAdicercall()] function and 
#' stored in the column `DicerConsensus` when \code{data.type="consensus"}. 
#' When plotting samples individually, there is the option to overlap the 
#' results onto a single line graph when \code{overlap=TRUE}. This is not an
#'option for bar plots.
#'
#'
#'
#'@return The function returns a list containing the results: a data frame and 
#'the plot(s). To access an element, simply use the "$" symbol,
#'and the elements "data" and "plot" will appear. 
#'
#' @examples
#' # load data 
#' data('sRNA_data')
#'
#' p1 <- RNAdistribution(data = sRNA_data, style = "line")
#'
#' p2 <- RNAdistribution(data = sRNA_data, style = "line", overlap = FALSE)
#'
#' p3 <- RNAdistribution(data = sRNA_data, style = "bar")
#'
#' p3.2 <- RNAdistribution(data = sRNA_data, style = "bar",
#'                        samples = c("heterograft_1", "heterograft_2",
#'                        "heterograft_3"))
#'
#' p4 <- RNAdistribution(data = sRNA_data, style = "bar", facet = FALSE)
#'
#' p5 <- RNAdistribution(data = sRNA_data, style = "bar",
#'                       facet = TRUE, facet.arrange = 2 )
#'
#'
#' # Run function to define sRNA class for each cluster.
#' sRNA_data_dicercall <- mobileRNA::RNAdicercall(data = sRNA_data, tidy=TRUE)
#'
#'p6 <- RNAdistribution(data = sRNA_data_dicercall, style = "bar", data.type = "consensus")
#'
#' @export
#' @importFrom BiocGenerics grep
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr count
#' @importFrom dplyr mutate
#' @importFrom tidyselect all_of
#' @importFrom data.table setDT
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 unit
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 margin
#' @importFrom tidyr gather
#' @importFrom data.table melt
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 expansion
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 scale_fill_manual
#'
RNAdistribution  <- function (data, samples = NULL, style, 
                              data.type = "samples", 
                              facet = TRUE, 
                              facet.arrange = 3, 
                              colour = "#0868AC", 
                              outline = "black", 
                              wrap.scales = "fixed", 
                              overlap = TRUE, 
                              relative = FALSE,
                              non.classified = TRUE) 
{
  if (base::missing(data) || !base::inherits(data, c("data.frame"))) {
    stop("data must be a data frame, see ?help for more details")
  }
  if (missing(style) || !is.character(style)) {
    stop("style parameter is missing or not a character vector.")
  }
  allowed_styles <- c("line", "bar")
  if (!style %in% allowed_styles) {
    stop("style parameter must be one of 'line', 'bar', or 'consensus'.")
  }
  allowed_datatypes <- c("samples", "consensus")
  if (!data.type %in% allowed_datatypes) {
    stop("data.type parameter must be one of 'samples' or 'consensus'.")
  }
  custom_theme <- ggplot2::theme(
    plot.title = ggplot2::element_text(face = "italic", size = 17), 
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(color="black", size = 15,
                                        face = "bold", 
                                        margin = ggplot2::margin(t = 10, b = 4)),
    axis.text.y = ggplot2::element_text(color="black", size = 15,face = "bold", 
                                        margin = ggplot2::margin(r = 10)) ,
    panel.grid.major.x = ggplot2::element_line( size=.1, color="grey", 
                                                linetype = 2 ),
    panel.grid.major.y = ggplot2::element_line( size=.1, color="grey", linetype = 2 ),
    legend.position = "right", 
    legend.box.margin= ggplot2::margin(20,20,20,20),
    legend.text = ggplot2::element_text(size=14, margin = ggplot2::margin(7,7,7,7)),
    legend.title = ggplot2::element_text(size = 14.5, face = "bold"),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10), size = 17, face = "bold"),
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10), size = 17, face = "bold"),
    strip.background = ggplot2::element_rect(fill = "lightgrey", colour = "black"), 
    strip.placement = "outside",
    strip.text.x = ggplot2::element_text(size = 14,face="italic" , margin = margin(b = 6, t = 5)),
    panel.border = ggplot2::element_rect(fill = "transparent", color = "black", linewidth = 1.5),
    plot.margin = ggplot2::unit(c(0.1, 0.1, 0.1, 0.1), "inches"))
  
  if (data.type == "consensus") {
    x <- data %>% dplyr::count(DicerConsensus)
    if (!relative == FALSE) {
      x <- x %>% dplyr::mutate(freq = n/sum(n))
      if (style == "line") {
        p1 <- ggplot2::ggplot(x, ggplot2::aes(x = DicerConsensus, 
                                              y = freq, group = 1)) + 
          ggplot2::geom_point(color = colour) + 
          ggplot2::geom_line(color = colour) + 
          ggplot2::theme_classic() + 
          ggplot2::xlab("sRNA Class") + ggplot2::ylab("Relative frequency") + 
          ggplot2::scale_y_continuous(expand = c(0, 0), 
                                      limits = c(0, 1)) + 
          ggplot2::theme_classic() + 
          custom_theme
        
      }
      else if (style == "bar") {
        
        p1 <- ggplot2::ggplot(x, ggplot2::aes(x = DicerConsensus, 
                                              y = freq, group = 1)) + 
          ggplot2::geom_bar(stat = "identity", fill = colour, color = outline) + 
          ggplot2::xlab("sRNA Class") + 
          ggplot2::ylab("Relative frequency") + 
          ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0,  0)) + 
          ggplot2::theme_classic() + 
          custom_theme
      }
    }
    else if (relative == FALSE) {
      perc_lim <- max(x$n) * 1.10
      
      if (style == "line") {
        max_lim <- max(x$n)
        
        p1 <- ggplot2::ggplot(x, ggplot2::aes(x = DicerConsensus, 
                                              y = n, group = 1)) + 
          ggplot2::geom_point(color = colour) + 
          ggplot2::geom_line(color = colour) + 
          ggplot2::xlab("sRNA Class") + ggplot2::ylab("Count") + 
          ggplot2::scale_y_continuous(expand = c(0, 0),limits = c(0, perc_lim))+ 
          ggplot2::theme_classic() + 
          custom_theme
      }
      else if (style == "bar") {
        p1 <- ggplot2::ggplot(x, ggplot2::aes(x = DicerConsensus, 
                                              y = n, group = 1)) + 
          ggplot2::geom_bar(stat = "identity", fill = colour, color = outline) + 
          ggplot2::xlab("sRNA Class") + 
          ggplot2::ylab("Count") + 
          ggplot2::scale_y_continuous(limits = c(0, perc_lim), expand = c(0,  0)) + 
          ggplot2::theme_classic() + 
          custom_theme
      }
    }
    out <- list(plot = p1, data = x)
    return(out)
  }
  else {
    if (base::missing(style) || !style %in% c("bar", "line")) {
      stop("Please specify type of plot to produce (line or bar)")
    }
    data.cols <- data %>% dplyr::select(tidyselect::starts_with("DicerCall_"))
    counts.df <- apply(data.cols, MARGIN = 2, table)
    if (base::inherits(counts.df, c("list"))) {
      class_colnames <- colnames(data)[grep("DicerCall_", 
                                            colnames(data))]
      required_columns <- unique(unlist(data[class_colnames]))
      for (i in 1:length(counts.df) ) {
        table_i <- counts.df[[i]]
        if (length(names(table_i)) < length(required_columns)) {
          missing_columns <- setdiff(required_columns, 
                                     names(table_i))
          table_i[missing_columns] <- 0
          if (!identical(names(table_i), required_columns)) {
            table_i <- table_i[match(required_columns, 
                                     names(table_i))]
          }
          counts.df[[i]] <- table_i
        }
      }
      # reorder dimensions
      desired_order <- c("N", "21", "22", "24", "23")
      counts.df <- lapply(counts.df, reorder_table, order = desired_order)
      
      counts.df <- data.frame(t(do.call(rbind, counts.df)))
    } else if (!base::inherits(counts.df, c("list"))) {
      counts.df <- data.frame(counts.df)
    }
    colnames(counts.df) <- gsub("DicerCall_", "", colnames(counts.df))
    counts.df <- data.table::setDT(counts.df, keep.rownames = "Class")[]
    if(non.classified == FALSE){
      counts.df <- counts.df[!(counts.df$Class == "N"), ]
    }
    if (style == "bar") {
      if (!is.null(samples)) {
        counts_class <- counts.df %>% select(Class)
        counts.df <- counts.df %>% select(all_of(samples))
        counts.df <- cbind(counts_class, counts.df)
      }
      plist <- list()
      
      for (col in names(counts.df)[-grep("Class", names(counts.df))]) {
        xlim_v <- max(counts.df[[col]])*1.10
        
        p <- if (!is.null(colour) && length(colour) == (nrow(counts.df) -1) ) {
          ggplot2::ggplot(counts.df, ggplot2::aes(x = .data[["Class"]],  y = .data[[col]]))+
            ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = variable), color = outline)+
            ggplot2::scale_fill_manual(values = colour )+
            ggplot2::theme_classic() + 
            ggplot2::labs(title = paste0("Sample: ", col), 
                          x = "sRNA Class", y = "Count") + 
            ggplot2::theme_bw() + 
            custom_theme
        } else if (!is.null(colour) && length(colour) == 1) {
          ggplot2::ggplot(as.data.frame(counts.df), ggplot2::aes(x = .data[["Class"]],  y = .data[[col]]))+
            ggplot2::geom_bar(stat = "identity", fill = colour, color = outline)+
            ggplot2::theme_classic() + 
            ggplot2::labs(title = paste0("Sample: ", col), 
                          x = "sRNA Class", y = "Count") + 
            ggplot2::theme_bw() + 
            custom_theme
        } else {
          ggplot2::ggplot(counts.df, ggplot2::aes(x = .data[["Class"]],  y = .data[[col]]))+
            ggplot2::geom_bar(stat = "identity",  ggplot2::aes(fill = variable), 
                              color = outline)+
            ggplot2::theme_classic() + 
            ggplot2::labs(title = paste0("Sample: ", col), 
                          x = "sRNA Class", y = "Count") + 
            ggplot2::theme_bw() + 
            custom_theme
        }
        
        
        
        plist[[col]] <- p
      }
      plist <- lapply(plist, function(p) {
        p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 
                                                                             0.1)))
      })
      sn <- names(plist)
      if (facet == TRUE) {
        p <- ggplot2::ggplot(tidyr::gather(counts.df, 
                                           key, Count, -Class), ggplot2::aes(Class, Count)) + 
          ggplot2::geom_bar(stat = "identity", fill = colour, 
                            colour = outline) + ggplot2::facet_wrap(~key, 
                                                                    scales = wrap.scales, ncol = facet.arrange) + 
          ggplot2::labs(x = "sRNA Class", y = "Count") + 
          ggplot2::theme_classic()+
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.1)))+ 
          custom_theme
        
        out <- list(plot = p, data = counts.df)
      }
      else if (facet == FALSE) {
        out <- list(plot = plist, data = counts.df)
      }
      return(out)
    }
    if (style == "line") {
      if (overlap == TRUE) {
        if (is.null(samples)) {
          counts.df_melt <- data.table::melt(counts.df, 
                                             id.vars = "Class")
          p <- ggplot2::ggplot(counts.df_melt, ggplot2::aes(Class, 
                                                            value, col = variable, group = 1)) + ggplot2::geom_point() + 
            ggplot2::geom_line() + ggplot2::xlab("sRNA Class") + 
            ggplot2::ylab("Count") + ggplot2::labs(color = "Samples") + 
            ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 21,size = 8))) + 
            ggplot2::theme_bw() + 
            custom_theme
          out <- list(plot = p, data = counts.df)
        }
        else if (!is.null(samples)) {
          counts.df <- counts.df %>% select(Class, all_of(samples))
          counts.df <- data.table::melt(counts.df, id.vars = "Class")
          p <- ggplot2::ggplot(counts.df, 
                               ggplot2::aes(Class, value, col = variable, group = 1)) + 
            ggplot2::geom_point(ggplot2::aes(color = variable)) + 
            ggplot2::geom_line(ggplot2::aes(color = variable)) +
            {if(!is.null(colour) && length(colour)== length(unique(counts.df$variable)))
              ggplot2::scale_colour_manual(values = colour ) }+
            ggplot2::xlab("sRNA Class") + 
            ggplot2::ylab("Count") + ggplot2::labs(color = "Samples") + 
            ggplot2::theme_bw() + 
            custom_theme
          
          out <- list(plot = p, data = counts.df)
        }
      }
      else if (overlap == FALSE) {
        if (facet == TRUE) {
          if (is.null(samples)) {
            p <- ggplot2::ggplot(tidyr::gather(counts.df, key, Count, -Class), 
                                 ggplot2::aes(Class, Count, group = 1)) + 
              ggplot2::geom_point(colour = colour) + 
              ggplot2::geom_line(colour = colour) + 
              ggplot2::facet_wrap(~key, scales = wrap.scales, 
                                  ncol = facet.arrange) + 
              ggplot2::xlab("sRNA Class") + 
              ggplot2::ylab("Count") + 
              ggplot2::theme_bw() + 
              ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) + 
              custom_theme                                                                                                                                                                                                                                                                                                                                            
          }
          else if (!is.null(samples)) {
            counts_class <- counts.df %>% select(Class)
            counts.df <- counts.df %>% select(all_of(samples))
            counts.df <- cbind(counts_class, counts.df)
            p <- ggplot2::ggplot(tidyr::gather(counts.df, 
                                               key, Count, -Class), ggplot2::aes(Class, 
                                                                                 Count, group = 1)) + ggplot2::geom_point(colour = colour) + 
              ggplot2::geom_line(colour = colour) + ggplot2::theme_classic() + 
              ggplot2::facet_wrap(~key, scales = wrap.scales, 
                                  ncol = facet.arrange) + ggplot2::theme_bw() + 
              ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.1))) +   
              custom_theme
          }
          out <- list(plot = p, data = counts.df)
          return(out)
        }
        else if (facet == FALSE) {
          if (!is.null(samples)) {
            counts_class <- counts.df %>% select(Class)
            counts.df <- counts.df %>% select(all_of(samples))
            counts.df <- cbind(counts_class, counts.df)
          }
          plist2 <- list()
          for (col in names(counts.df)[-grep("Class",  names(counts.df))]) {
            
            xlim_v <- max(counts.df[[col]])*1.10
            
            p2 <- if (!is.null(colour) && length(colour) == (nrow(counts.df) -1) ) {
              ggplot2::ggplot(counts.df, 
                              ggplot2::aes(x = .data[["Class"]],  y = .data[[col]], group = 1))+
                ggplot2::geom_point(ggplot2::aes(color = variable))+
                ggplot2::geom_line(ggplot2::aes(color = variable)) + 
                ggplot2::scale_fill_manual(values = colour )+
                ggplot2::theme_classic() + 
                ggplot2::labs(title = paste0("Sample: ", col), 
                              x = "sRNA Class", y = "Count") + 
                ggplot2::theme_bw() + 
                custom_theme
            } else if (!is.null(colour) && length(colour) == 1) {
              ggplot2::ggplot(as.data.frame(counts.df), 
                              ggplot2::aes(x = .data[["Class"]],  y = .data[[col]], group = 1))+
                ggplot2::geom_point( color = colour)+
                ggplot2::geom_line(color = colour) + 
                ggplot2::theme_classic() + 
                ggplot2::labs(title = paste0("Sample: ", col), 
                              x = "sRNA Class", y = "Count") + 
                ggplot2::theme_bw() + 
                custom_theme
            } else {
              ggplot2::ggplot(counts.df, 
                              ggplot2::aes(x = .data[["Class"]],  y = .data[[col]], group = 1))+
                ggplot2::geom_point(ggplot2::aes(color = variable))+
                ggplot2::geom_line(ggplot2::aes(color = variable)) + 
                ggplot2::theme_classic() + 
                ggplot2::labs(title = paste0("Sample: ", col), 
                              x = "sRNA Class", y = "Count") + 
                ggplot2::theme_bw() + 
                custom_theme
            }
            
            plist2[[col]] <- p2
          }
          plist2 <- lapply(plist2, function(p2) {
            p2 + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 
                                                                                  0.1)))
          })
          sn <- names(plist2)
          p <- plist2
          out <- list(plot = p, data = counts.df)
        }
      }
      return(out)
    }
  }
}