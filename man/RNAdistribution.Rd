% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNAdistribution.R
\name{RNAdistribution}
\alias{RNAdistribution}
\title{Plot the distribution of sRNA lengths}
\usage{
RNAdistribution(
  data,
  samples = NULL,
  style,
  facet = TRUE,
  facet.arrange = 3,
  colour = "black",
  together = TRUE,
  consensus = FALSE,
  relative = FALSE
)
}
\arguments{
\item{data}{a dataframe, on which one of the following functions has already
been called: \code{\link[=RNAimport]{RNAimport()}},\code{\link[=RNAdicercall]{RNAdicercall()}}.}

\item{samples}{character vector. Store names of samples to analyse and plot.
Argument is required for plotting individual sample replicates, either
individually, overlapped together or in a facet. Use the sample
replicate names present in the data frame (\code{data}), select samples you wish
to plot. Is not required when plotting the sRNA consensus using the argument
\code{total=TRUE}.}

\item{style}{plotting option to choose the style of either a line graph or
bar chart to represent your data.
\itemize{
\item Where \code{style="line"} a line graph will be produced
\item Where \code{style="bar"} produces a bar graph
\item Where \code{style="consensus"} produces the line graph for the consensus
sRNA in conjunction with \code{consensus=TRUE}
}}

\item{facet}{Logical; forms a matrix of panels defined by row and column
faceting variables. It plots the results for each sample as a bar-chart
and contains it within a single plot. The number of rows in the facet can
be changed using the argument \code{facet.arrange}.
Default \code{facet = TRUE} , plots each sample separately when
\code{facet = FALSE} .}

\item{facet.arrange}{numeric; value supplied to define the number of columns
to include in the facet. This argument is piped into the \code{ ncol}  argument in
\code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}} to define the number of columns. By default, this is
set to 3.}

\item{colour}{bar plot fill colour. Default colour is "black".}

\item{together}{Logical; forms a single line graph with multiple lines each
to represent the sample replicates. Default \code{together=TRUE}.}

\item{consensus}{Logical; plots the distribution of sRNA classes across all
identified dicer-derived clusters based on the consensus. See
\code{\link[=RNAdicercall]{RNAdicercall()}}
function to calculate consensus RNA class based the experimental replicates.
Default \code{consensus=FALSE}.}

\item{relative}{Logical; calculates relative frequency of consensus
dicer-derived sRNA clusters. Only applicable when  only in conjunction with
argument \code{consensus}, when set as \code{consensus=TRUE}. As default,
\code{relative=FALSE}.}
}
\value{
The function returns a list containing the results in the form of a
data frame and the plot(s). To access an element, simply use the "$" symbol,
and the elements "data" and "plot" will appear. The \code{samples} argument allows
uses to plot specific samples in a single plot (facet bar plot or line graph).
This can encourage closer comparison between sample replicates.
}
\description{
\code{RNAdistribution} plots the distribution of dicer-derived
sRNA classes (20-24nt) across samples or across the sRNA consensus
determined by the function \code{\link[=RNAdicercall]{RNAdicercall()}}.
}
\details{
The function can be used to plot a variety of different comparisons and plots.
It can be used to plot the distribution of sRNA classes within each sample
replicate, which can be represented as a bar chart \code{style="bar"} or a
line graph \code{style="line"}. These plots can be represented individually or
in a single plot facet \code{facet="TRUE"} by default.

To plot the sRNA dicer-derived clusters identified in each sample, the
function extracts the information from the RNA summary data and calculates the
total number of each RNA class identified within a sample, for all samples.

Alternatively, the function allows you to plot the line graph for each sample
together, overlapped on a single graph \code{total="TRUE"}. This is not an
option for bar plots.

The final option, is to plot the total consensus of dicer-derived sRNA
clusters across the experimental conditions, the function pulls the consensus
call from the column created by the \code{\link[=RNAdicercall]{RNAdicercall()}} function in the working
data frame.
}
\examples{
data('sRNA_data')

p1 <- RNAdistribution(data = sRNA_data, style = "line")

p1.2 <- RNAdistribution(data = sRNA_data, style = "line",
                        samples = c("heterograft_1", "heterograft_2",
                        "heterograft_3"))
p2 <- RNAdistribution(data = sRNA_data, style = "line", together =FALSE )

p3 <- RNAdistribution(data = sRNA_data, style = "bar")

p3.2 <- RNAdistribution(data = sRNA_data, style = "bar",
                       samples = c("heterograft_1", "heterograft_2",
                       "heterograft_3"))

p4 <- RNAdistribution(data = sRNA_data, style = "bar", facet = FALSE)

p5 <- RNAdistribution(data = sRNA_data, style = "bar",
                      facet = TRUE, facet.arrange = 2 )

data("sRNA_data_consensus")
p6 <- RNAdistribution(data = sRNA_data_consensus, consensus = TRUE)

}
