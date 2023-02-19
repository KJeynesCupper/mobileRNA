#' Calculate chromosome sizes, and store in a text file.
#'
#' @description Invokes an OS command to calculate the chromosome sizes within
#' a genome.
#' @param  genome Directory path to the genome (\code{.fa/.fasta})
#' @param  out Path to output directory
#' @return A plain text file containing the chromosome sizes
#' @examples
#' \dontrun{
#' ref <- "/Users/user1/projectname/workplace/reference/cat/ref_merged.fa"
#' save <- "/Users/user1/projectname/workplace/reference/cat/"
#' chr_sizes(ref)
#' }
#'
chr_info <- function(genome, out){
  system(paste0("samtools faidx ", genome ))
  system(paste0( "cut -f1,2 ", genome,".fai > ", out,"chr_sizes.txt"))
}
