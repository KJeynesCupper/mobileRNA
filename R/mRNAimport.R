#' Import and organise mRNAseq data set
#'
#' @description Organises mRNAseq counts data from each sample replicate into
#' a single dataframe.
#'
#'
#' @details
#' Supply the directory to the folder storing the plain text files for each
#' sample replicate contain counts data. Each file is draw from the directory
#' based on the file names supplied to the samples argument, minus the .txt
#' extension. The file names supplied will represent the sample across the
#' analysis, unless changed or by altering column names after importation with
#' this function.
#'
#' The counts data is composed of two columns; the first containing gene names,
#' and the second containing counts.
#'
#' The first column of the function output contains all the gene names identified
#' across all samples. The other columns represent the counts for each sample
#' replicate imported. Where the sample name is incorporated into the column name.
#' For example, for a sample called "Sample1", the counts will be stored in a
#' column called "Count_Sample1".
#'
#' @param loci dataframe; rows represent individual sRNA dicer-derived clusters
#' columns state genomic location (labeled "Locus") and the cluster name
#' (labeled "Cluster"). This file contains an accumulation of all
#' sRNA dicer-derived clusters  identified across samples in the analyse, without
#' duplication. This can be produced using the [mobileRNA::RNAloci()]
#' function
#'
#'
#'
#' @param directory Path to directory containing of sample folders.
#'
#' @param samples Vector of characters naming the sample names correlating
#' to outputted folders located in the \code{`directory`} argument path.
#'
#'
#' @param report Logical; prompts progress.
#'
#'@param tidy Logical; removes genes from analysis where there are zero counts
#'across all replicates.
#'
#'
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @importFrom data.table fread
#' @importFrom data.table setnames
#' @importFrom dplyr mutate across contains
#' @importFrom tidyr replace_na
#' @importFrom magrittr %>% filter
#' @importFrom stats unique
#' @importFrom utils cat
#' @importFrom methods mget
mRNAimport <- function(directory, samples, report = TRUE, tidy = TRUE){

  # load data as list
  sample_data <- list()
  for (file in samples) {
    sample_data[[file]] <- data.table::fread(paste0(directory, file, ".txt"), header = FALSE)
    colnames(sample_data[[file]])[1] <- "Gene"
    colnames(sample_data[[file]])[2] <- "Count"
  }

  # check each file has two columns
  for (df in sample_data) {
    if (!all(ncol(df) == 2)) {
      stop("mRNA dataset(s) does not contain required columns")
    }
  }
  # merge first columns to create list of genes across all samples
  genes <- lapply(sample_data, "[", , "Gene")
  genes_all <- unique(Reduce(merge,genes))

  # ADDs sample information to the genes_all object
  for (i in seq_along(sample_data)){
    matches <- genes_all[sample_data[[i]], on = "Gene", nomatch = 0]
    matches_values <- matches[, .(Count=sum(Count)),by = "Gene"]
    # Rename the aggregated columns
    col_name <- paste0("Count_", names(sample_data)[i])
    data.table::setnames(matches_values, c("Gene", col_name))
    # Merge the aggregated values back into df1
    genes_all[matches_values, on = "Gene", (col_name) := mget(col_name)]
    # Print progress
    if (report) {
      cat("Added information from", i, "to the analysis mRNA dataframe" ,"\n")
    }
  }


  # Fill in missing values with 0 or N
  mRNA_information <- genes_all %>%
    dplyr::mutate(dplyr::across(dplyr::contains('Count_'), ~tidyr::replace_na(.,0)))
  # Convert data.frame and return it
  mRNA_information <- as.data.frame(mRNA_information)

  # remove rows with
  if (tidy){
    mRNA_information <- mRNA_information %>%
      filter(if_any(where(is.numeric), ~. != 0))
    }

  return(mRNA_information)
}

