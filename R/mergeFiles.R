#' Merge two reference genomes into one
#'
#' @description Merges two reference genomes (`.fa/.fasta)`. into one single
#' reference.
#'
#'
#' @param  files Path to the directory containing genome both reference files
#' (The directory must only contain these files).
#' @param  file1 Path to the directory of the first genome reference file.
#' @param  file2 Path to the directory of the second genome reference file.
#' @param  out Path desired output directory.
#'
#'
#' @return A file named \code{ref_merged.fa} representing the merged
#' genome reference.
#'
#' @details Prior to merging it is import to check the chromosome names for each
#' genome. This is to ensure they can be easily distinguished from each other
#' after merging. To extract the chromosomal information, look at
#' [RNAlocate::referenceInfo()] and to modify the chromosome names, look at
#' [RNAlocate::modifyReference()].
#'
#'
#' @examples \dontrun{
#' mergeFiles(file1 = "./workplace/reference/ref1.fa.gz",
#'             file2 = "./workplace/reference/ref2.fa.gz"
#'             out = "./workplace/reference/merge/ref_merged.fa")
#' # Or,
#' mergeFiles(files = "./workplace/reference/*",
#'             out = "./workplace/reference/merge/ref_merged.fa")
#'}
#'
#'
#' @importFrom stringr "str_detect"
#' @export
mergeFiles <- function(files, file1=NULL, file2=NULL, out){
  system(paste0("gunzip " ,files,"*.gz"))
  if(!is.null(file1) | !is.null(file2)) {
    if(stringr::str_detect(file1, ".fa") || stringr::str_detect(file1,".fasta")
       & stringr::str_detect(file2, ".fa") ||
       stringr::str_detect(file2, ".fasta")) {
      system(paste0("gunzip " ,file1))
      system(paste0("gunzip " ,file2))
      system(paste0("cat ", file1, " ", file2, " > ", out))
    } else
      if(stringr::str_detect(file1, ".gff") || stringr::str_detect(file1,
                                                                   ".gff3")
         & stringr::str_detect(file2, ".gff") ||
         stringr::str_detect(file2, ".gff3")) {
        system(paste0("gunzip " ,file1))
        system(paste0("gunzip " ,file2))
        system(paste0("cat ", file1, " ", file2, " > ", out))
  }
    if(stringr::str_detect(files, ".fa") ||
       stringr::str_detect(files, ".fasta")) {
    system(paste0("cat ", files," > ", out))
    }
    if(stringr::str_detect(files, ".gff") ||
       stringr::str_detect(files, ".gff3")) {
      system(paste0("cat ", files," > ", out))
    }
  }
}




