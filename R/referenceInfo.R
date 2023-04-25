#' Extract chromosome sizes
#'
#' @description Invokes an OS command to calculate the chromosome sizes within
#' a given genome.
#'
#'
#' @param  genome Path to a genome (\code{.fa/.fasta})
#'
#'
#' @param out Path to output directory
#'
#' @return A plain text file listing the chromosomes and their respective sizes.
#'
#'
#' @examples
#' \dontrun{
#' ref <- "/Users/user1/projectname/workplace/reference/merge/ref_merged.fa"
#' save <- "/Users/user1/projectname/workplace/reference/merge/"
#' chrInfo(genome = ref, out = save)
#' }
#'
referenceInfo <- function(genome, out){
  system(paste0("samtools faidx ", genome ))
  system(paste0( "cut -f1,2 ", genome,".fai > ", out,"chr_sizes.txt"))
}


