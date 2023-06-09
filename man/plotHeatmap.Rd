% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotHeatmap.R
\name{plotHeatmap}
\alias{plotHeatmap}
\title{Heatmap using hierarchical clustering}
\usage{
plotHeatmap(data, colours = NULL, dendogram = TRUE, cellheight = NULL)
}
\arguments{
\item{data}{Data frame containing FPKM values for each samples on
each sRNA dicer-derived cluster of interest.}

\item{colours}{Colours to display and represent the heatmap.
Defaults to \code{\link[grDevices:palettes]{grDevices::heat.colors()}} (heat.colors(100)).}

\item{dendogram}{Logical; indicating whether to include the dendrogram and
clustering, and retain clustering. Default, \code{dendogram = TRUE} to
include.}

\item{cellheight}{individual cell height in points. If left as NA, then the
values depend on the size of plotting window.}
}
\value{
Produces a list of objects, including the plot.
}
\description{
Plots a heatmap with hierarchical clustering via an rlog
transformation of RPM data and euclidean statistics.
}
\details{
The function create a heatmap based on the hierarchical clustering
of FPKM values using euclidean statistics.
}
\examples{


data("sRNA_data_mobile")

# plot heatmap of potential mobile sRNAs
 p1 <-  plotHeatmap(sRNA_data_mobile)




}
