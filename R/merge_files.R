#' Merge the reference genomes into a single reference genome
#'
#' @description Merges two reference genomes (`.fa/.fasta)`. into one.
#' @param  files Path to the directory containing genome both reference files
#' (The directory must only contain these files).
#' @param  file1 Path to the directory of the first genome reference file.
#' @param  file2 Path to the directory of the second genome reference file.
#' @param  out Path desired output directory.
#' @return A file named `"ref_merged.fa"` representing the merged
#' genome reference.
#' @examples \dontrun{
#' merge_files(file1 = "./workplace/reference/ref1.fa.gz",
#'             file2 = "./workplace/reference/ref2.fa.gz"
#'             out = "./workplace/reference/merge/")
#' # Or,
#' merge_files(files = "./workplace/reference/",
#'             out = "./workplace/reference/merge/")
#'}
#'
#' @export
merge_files <- function(files, file1=NULL, file2=NULL, out){
  system(paste0("gunzip " ,files,"*.gz"))
  if(!is.null(file1) | !is.null(file2)) {
    system(paste0("gunzip " ,file1))
    system(paste0("gunzip " ,file2))
    system(paste0("cat ", file1, " ", file2, " > ", out, "ref_merged.fa"))
  } else {
    system(paste0("cat ", files," > ", out, "ref_merged.fa"))
  }
}


