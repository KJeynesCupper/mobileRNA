#' Remove punctuation from chromosome names.
#'
#' @description Invokes an OS command to alter the names of chromosomes by
#' selecting a character sting within the chromosomes names to replace.
#' See [funs(chr_checker)] to identify the names of the chromosomes in a genome.
#' @param  old Character string of chromosome name to change.
#' @param  new Replacement of old string. A character string to rename the
#' chromosomes.
#' @param  genome Directory to the location of the merged genome reference file.
#' @return Edited ref_merged.fasta file representing the merged genome
#' reference, saved in the "cat" directory within "reference" directory.
#' @examples
#' \dontrun{
#' # The Solanum lycopersicum (version 4) genome contains a
#' # full-stop/period within each chromosome name which needs to be removed.
#' # Each chromosome name starts with "SL4." and we will be changing it to
#' # remove the full-stop.
#'
#'
#' old <- "SL4."
#' new <- "SL4"
#' genome <- "/Users/user1/projectname/workplace/reference/cat/ref_merged.fa"
#' chrModify(old, new, genome)
#' }
#'
#' @export
chrModify <- function(old, new, genome){
  system(paste0("gunzip " , genome ))
  system(paste0("sed -i ","'s/",old,"/",new,"/g"," ", genome))
}

