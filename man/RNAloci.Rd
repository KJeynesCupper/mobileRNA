% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNAloci.R
\name{RNAloci}
\alias{RNAloci}
\title{Compute complete list of sRNA dicer-derived cluster loci}
\usage{
RNAloci(files, out, samples)
}
\arguments{
\item{files}{Path to the directory containing ShortStack mapping results.}

\item{out}{Path to the directory to save the outputted annotation files,
include the name of the output file, including file extension (.txt`)}

\item{samples}{A character vector containing the names of all the samples in
the analysis. These should match the sample folder names created in mapping
step one, outputted by \code{Shortstack}.}
}
\value{
A plain text file (.txt) and when assigned to a variable, the
dataframe is saved to the global environment. The output contains two
columns,\code{Locus} containing the genomic coordinates of the dicer-derived
cluster and \code{Cluster} contains the given name for the dicer-derived cluster.
}
\description{
Collates all in the identified sRNA dicer-derived clusters
across multiple sample replicates.
}
\details{
Based on the output of Shortstack, the function pulls the loci
of the dicer-derived clusters within each sample, supplied to the function,
and merges the loci into a single data frame.

The function utilises the \code{ShortStack_All.gff3} file produced from
the mapping and clustering analysis by ShortStack.

The function outputs the data frame to the global environment when assigned
to a variable/object, but regardless of this it outputs the data frame as a
plain text file in the given directory.
}
\examples{
\dontrun{


folder <- "./workplace/alignment_unique_one/"
save_folder <- "./workplace/reference/LociInfo.txt"

conditions <- c("heterograft_1","heterograft_2","heterograft_3",
                 "selfgraft_1","selfgraft_2", "selfgraft_3")


Loci_info <- RNAloci(files = folder,
                    out = save_folder,
                    samples = conditions)

}
}
