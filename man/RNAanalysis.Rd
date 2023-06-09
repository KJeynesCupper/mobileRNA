% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNAanalysis.R
\name{RNAanalysis}
\alias{RNAanalysis}
\title{Differential Expression (DE) Analysis of sRNA dicer-derived clusters
using \code{DESeq2} or \code{edgeR}}
\usage{
RNAanalysis(data, group, method = c("edgeR", "DESeq2"), dispersionValue = NULL)
}
\arguments{
\item{data}{numeric data frame produced by \code{\link[=RNAimport]{RNAimport()}} and/or
\code{\link[=RNAsubset]{RNAsubset()}}.}

\item{group}{Vector of the condition (ie. treatment or control) for each
sample. Must be stated in the same order as the samples in the \code{data} file
from left to right.}

\item{method}{The method to undertaken differential analysis, choose from
methods of either \link[DESeq2:DESeq]{DESeq2::DESeq} or \link[edgeR:edgeR-package]{edgeR::edgeR}.
Must be stated as either "DESeq2" or "edgeR" in the function.}

\item{dispersionValue}{numeric value; manual setting of dispersion value
which is recommended for analysis in experiments without biological
replicates when utilising the \link[edgeR:edgeR-package]{edgeR::edgeR} method.}
}
\value{
Undertakes differential analysis, based on a specified method, and
adds the results to the supplied dataframe. This includes:
\itemize{
\item Count mean
\item Log fold change
\item p-value
\item Adjusted p-value
\item Log counts per million (CPM/RPM)
}
}
\description{
This function allows you to compute the differential expression
of sRNA dicer-derived clusters. The function allows the choice between
analysis with \code{DESeq2} or \code{edgeR}.
}
\details{
The analysis allows the users to choose the method which best suits
their data. Notably, \code{DESeq2} cannot compute the analysis when there only
one replicate per condition, but, \code{edgeR} can. Simply set a suitable
dispersion value, based on similar data, to use this feature. The dispersion
value is other wise known as the common Biological squared coefficient
of variation. A typical dispersion value is 0.4 for human data sets, 0.1 for
data on genetically identical model organisms or 0.01 for technical
replicate. See the User’s Guide for the  ‘EdgeR’ package for more details,
\link[edgeR:edgeR-package]{edgeR::edgeR}.
}
\examples{

data("sRNA_data_consensus")
# sample conditions.
groups <- c("Heterograft", "Heterograft", "Heterograft",
          "Selfgraft", "Selfgraft", "Selfgraft")


## Differential analysis: DEseq2 method
sRNA_DESeq2 <- RNAanalysis(data = sRNA_data_consensus,
                             group = groups,
                             method = "DESeq2" )


## Differential analysis: edgeR method
sRNA_edgeR <- RNAanalysis(data = sRNA_data_consensus,
                            group = groups,
                            method = "edgeR" )
}
