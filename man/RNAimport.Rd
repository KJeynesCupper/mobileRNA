% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNAimport.R
\name{RNAimport}
\alias{RNAimport}
\title{Import and organise sRNAseq or mRNA data sets}
\usage{
RNAimport(input = c("sRNA", "mRNA"), directory, samples, tidy = TRUE)
}
\arguments{
\item{input}{string; define type of Next-Generation Sequencing dataset
originates from, either "sRNA" or "mRNA" are the only valid inputs.}

\item{directory}{Path to directory containing of sample folders. NOTE:
Following the suggested pre-processing steps, these can be found in second
alignment folders.}

\item{samples}{Vector of characters naming the sample names correlating
to outputted folders located in the \code{`directory`} argument path.}

\item{tidy}{Logical; removes genes from analysis where there are zero counts
across all replicates.}
}
\value{
A dataframe where rows represent sRNA clusters and columns represent
replicate information. Replicate information includes Dicercall, Counts, RPM
and MajorRNA sequence. Each replicate information is distinguishable as
the replicate name is ajoined as a prefix to each column name.
}
\description{
Load and organise either sRNAseq or mRNAseq results into a
single dataframe containing all experimental replicates specified where rows
represent either a sRNA locus or gene, respectively.
}
\details{
Supply the directory to the folder storing the various sample replicate
folders produced by ShortStack analysis. Within each sample replicate folder,
there is "Result.txt". This file contains many different columns, but
for the analysis the columns of interest are  "Locus", "DicerCall", "Reads",
and "RPM" which will all be imported and loads into the new data frame.
Note that "Reads" will be converted to "Counts" as output of the function.
Locus contains the genomic locations of clusters,
Dicercall contains the most likely size of the cluster, Reads contains the
number of reads mapped to the cluster (recommended to use unique mapping,
hence these will be the number of uniquely mapped reads, ie, not including
multimapped reads). Lastly, RPM contains the Reads per Million score.

\code{"RNAimport"} takes the results from  pre-processed sRNAseq data which
has been mapped and undertaken cluster analysis through ShortStack, organises
into a single data frame. Each row represent an sRNA dicer-derived cluster
within the analysis and columns 1 - 5 supply information on the cluster
including the locus, the separated coordinates
(chromosome, start, end) and cluster name.

Further columns represent data imported for each samples including DicerCall,
Counts and RPM. The DicerCall represents the size of most
abundant small RNA size based on the parameter used in ShortStack.
The Count column represents the number of aligned sRNA-seq reads that overlap
the locus. The RPM represents the reads per million.
For each replicate included in the analysis, these columns are labeled with
the type and then then name of the sample, for example, for
a sample called "Sample1", the information from this sample will be stored in
columns DicerCall_Sample1, Count_Sample1 and RPM_Sample1.
}
\examples{
\dontrun{

# import sRNAseq data
df_sRNA <- RNAimport(input = "sRNA",
                     directory = "./analysis/sRNA_mapping_results/",
                     samples = c("heterograft_1", "heterograft_2",
                     "heterograft_3","selfgraft_1" , "selfgraft_2" ,
                     "selfgraft_3"))

# import mRNAseq data
df_mRNA <- RNAimport(input = "mRNA",
                     directory = "./analysis/mRNA_counts/",
                     samples = c("heterograft_1", "heterograft_2",
                     "heterograft_3","selfgraft_1" , "selfgraft_2" ,
                     "selfgraft_3"))

}
# The output of this function can be explore in the data object sRNA_data
data("sRNA_data")
head(sRNA_data)


}
