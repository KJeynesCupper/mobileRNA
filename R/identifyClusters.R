#' Merge the dicer-derived information into a single annotation file
#'
#' @description Collate all in the identified dicer-derived clusters, which
#' define each small RNA molecules in the analysis, into a single file.
#'
#' @param files The path to the local directory containing the results of
#' mapping step one.
#' @param out The path to the local directory to save the annotation files.
#' @param samples A character vector containing the names of all the samples in
#' the analysis. These should match the sample folder names.
#'
#' @return A plain text files (.txt) and General Feature Format (.gff3) file
#' containing the names of all cluster identified in across the samples
#' and the genomic locations.
#' @examples
#' \dontrun{
#' conditions <- c("TomEgg_1","TomEgg_2","TomEgg_3",
#'                  "TomTom_1","TomTom_2", "TomTom_3")
#'
#'folder <- "./workplace/alignment_multi_one/"
#'save_folder <- "./workplace/reference/"
#'
#'identifyClusters(files = folder,
#'            out = save_folder,
#'             samples = conditions)
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

identifyClusters <- function(files, out, samples ){
  gff_alignment <- GenomicRanges::GRangesList()
  for (i in samples){
    gff_alignment[[i]] <- import.gff(paste0(files,i, "/ShortStack_All.gff3"))
  }

  gff_merged <- GenomicRanges:: reduce(unlist(gff_alignment), ignore.strand=TRUE)

  gff_merged_df <- data.frame(position = paste0(as.character(
    GenomeInfoDb::seqnames(gff_merged)),":",stats::start(gff_merged),"-",
    stats::end(gff_merged)),name = paste0("cluster_", 1:length(gff_merged)))

  utils::write.table(gff_merged_df, file = paste0(out, "ClustersInfo.txt"),
              quote = F, sep = "\t", row.names = F, col.names = F)

  export(gff_merged, paste0(out, "ClustersInfo.txt"))

}
