% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNAsignificant.R
\name{RNAsignificant}
\alias{RNAsignificant}
\title{Extract statistically significant sRNA clusters}
\usage{
RNAsignificant(data, statistical = FALSE, padj = 0.05, p.value = NULL)
}
\arguments{
\item{data}{Numeric data frame}

\item{statistical}{If TRUE, will undertake statistical filtering based on the
a p-value or adjusted p-value threshold stated by \code{padj} & \code{p.value}.Default
set at FALSE. Requires presence of columns containing statistical data.
In order to filter by the adjusted p-value, a column named \code{padjusted} must
be present. Similarly, to filter by the p-value, a column named \code{pvalue} must
be present. See \code{\link[=RNAanalysis]{RNAanalysis()}} to calculate statistical values.}

\item{padj}{A user defined numeric value to represent the adjusted p-value
threshold to define statistic significance. Defaults set at 0.05.Only mobile
molecules with adjusted p-values equal or lower than specified are returned.}

\item{p.value}{A user defined numeric value to represent the p-value
threshold to define statistic significance. There is no default value, set
this instead of using an adjusted p-value to filter molecules. Only mobile
molecules with p-values equal or lower than specified are returned.}
}
\value{
A refined version of the working dataframe supplied to the function.
The function selects sRNA clusters which meet the statistical threshold,
given the statistical analysis has been undertaken using the
\code{\link[=RNAanalysis]{RNAanalysis()}} function.
}
\description{
Based on a given threshold, filter dataset for statistically
significant sRNA clusters.
}
\details{
Filters the data based on statistical significance
which has been calculated by the \code{\link[=RNAanalysis]{RNAanalysis()}} function. This
optional features enables the user to select sRNA clusters which meet a
specific p-value or adjusted p-values threshold.

Requires either or both columns: \code{pvalue}, \code{padjusted} to undertake.
}
\examples{

data("sRNA_data_consensus")

## sample conditions in order within dataframe
groups <- c("Heterograft", "Heterograft", "Heterograft",
            "Selfgraft", "Selfgraft", "Selfgraft")


## Differential analysis using the DEseq2 method 
sRNA_DESeq2 <- RNAanalysis(data = sRNA_data_consensus,
                          group = groups,
                          method = "DESeq2")
                          
 ## Select significant based on padjusted                        
significant_sRNAs <- RNAsignificant(sRNA_DESeq2)
                                 
                                 

}
