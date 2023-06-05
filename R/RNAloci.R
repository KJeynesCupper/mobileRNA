#' Compute complete list of sRNA dicer-derived cluster loci
#'
#' @description Collates all in the identified sRNA dicer-derived clusters
#' across multiple sample replicates.
#'
#' @details Based on the output of Shortstack, the function pulls the loci
#' of the dicer-derived clusters within each sample, supplied to the function, and
#' merges the loci into a single data frame.
#'
#' The function utilises the \code{ShortStack_All.gff3} file produced from
#' the mapping and clustering analysis by ShortStack.
#'
#' The function outputs the data frame to the global environment when assigned
#' to a variable/object, but regardless of this it outputs the data frame as a
#' plain text file in the given directory.
#'
#'
#'
#'
#'
#' @param files Path to the directory containing ShortStack mapping results.
#'
#' @param out  Path to the directory to save the outputted annotation files,
#' include the name of the output file, including file extension (.txt`)
#'
#' @param samples A character vector containing the names of all the samples in
#' the analysis. These should match the sample folder names created in mapping
#' step one, outputted by \code{Shortstack}.
#'
#' @return A plain text file (.txt) and when assigned to a variable, the dataframe
#' is saved to the global environment. The output contains two columns,
#' `Locus` containing the genomic coordinates of the dicer-derived cluster and
#' `Cluster` contains the given name for the dicer-derived cluster.
#'
#'
#' @examples
#' \dontrun{
#'
#'
#'folder <- "./workplace/alignment_unique_one/"
#'save_folder <- "./workplace/reference/LociInfo.txt"
#'
#'conditions <- c("TomEgg_1","TomEgg_2","TomEgg_3",
#'                  "TomTom_1","TomTom_2", "TomTom_3")
#'
#'
#'Loci_info <- RNAloci(files = folder,
#'                          out = save_folder,
#'                          samples = conditions)
#'
#' }
#' @export
#' @importFrom GenomicRanges "reduce"
#' @importFrom rtracklayer "import.gff"
#' @importFrom rtracklayer "export"
#' @importFrom GenomicRanges "GRangesList"
#' @importFrom GenomeInfoDb "seqnames"
#' @importFrom stats "start"
#' @importFrom stats "end"
#' @importFrom utils "write.table"

RNAloci <- function(files, out, samples ){

  if (base::missing(files)) {
    stop(paste("Please specify files, a connection to a local directory containing
               sample folders"))
  }
  if (base::missing(samples) || !base::inherits(samples, c("character"))) {
    stop(paste("Please specify samples, a vector containing individual strings
               cooresponding to the folders of each sample replicate containing
               results"))
  }
  if (base::missing(out) || tools::file_ext(out == "txt")) {
    stop(paste("Please specify out, a connection to a local directory to
               store output, including name and file extention (.txt)"))
  }

  gff_alignment <- GenomicRanges::GRangesList()
  for (i in samples){
    gff_alignment[[i]] <- rtracklayer::import.gff(paste0(files,i,
                                                      "/ShortStack_All.gff3"))
  }

  gff_merged <- GenomicRanges::reduce(unlist(gff_alignment), ignore.strand=TRUE)

  gff_merged_df <- data.frame(Locus = paste0(as.character(
    GenomeInfoDb::seqnames(gff_merged)),":",
    stats::start(gff_merged),
    "-",stats::end(gff_merged)),
    Cluster = paste0("cluster_", 1:length(gff_merged)))


  utils::write.table(gff_merged_df, file = out,
              quote = F, sep = "\t", row.names = F, col.names = T)
  message("Writting Loci file to:  ", out)
  return(gff_merged_df)
  message("Loci data frame saved to named object")
}
