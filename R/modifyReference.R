#' Modify chromosome names with a reference genome.
#'
#' @description Invokes an OS command to alter the names of chromosomes in a
#' given reference genome (.fa or .fasta) by selecting a character sting
#' within the chromosomes name to replace.
#'
#' See [RNAlocate::referenceInfo()] to identify the names of the chromosomes in a
#' genome.
#'
#'
#' @param  old Character string found across all chromosome in genome to change.
#' @param  new Replacement of old string. A character string to rename the
#' chromosomes.
#' @param  genome Directory to the location of the merged genome reference file.
#'
#' @details Strings should be contain within speech marks ("" or ''). This
#' function will run through all names in the reference and alter them.
#'
#' For the purpose of merging two genomes, it is important to make sure each
#' genome has distinguishable chromosome names from each other. In addition,
#' some clustering techniques require the chromosome names to be free of
#' punctuation. Therefore, it is recommended to remove any punctuation
#' from both reference genomes.
#'
#' When removing punctuation, ensure that other characters are included as it
#' will remove the punctuation universally. This could cause issues for
#' instance, if full-stops/periods were removed universally, as this could alter
#' numerical values.
#'
#' You can alter the chromosome names within reference genomes before or after
#' merging depending on the purpose. To make the names distinguishable from one
#' another, for example if they both use "NC" as prefixes, the chromosomes
#' should be modified before merging the genomes. Alternatively, if it is to
#' remove punctuation then this can be completed after merging genomes.
#'
#' @return Edits the existing genome reference file supplied.
#'
#'
#'
#' @examples
#' \dontrun{
#' # The Solanum lycopersicum (version 4) genome contains a
#' # full-stop/period within each chromosome name which needs to be removed.
#' # Each chromosome name starts with "SL4." and we will be changing it to
#' # remove the full-stop.
#'
#'
#' old_syn <- "SL4."
#' new_syn <- "SL4"
#' genome_ref <- "/Users/user1/projectname/workplace/reference/merge/ref_merged.fa"
#' chrModify(old = old_syn, new = new_syn, genome = genome_ref)
#' }
#'
#' @export
modifyReference <- function(old, new, genome){
  system(paste0("gunzip " , genome ))
  system(paste0("sed -i ","'s/",old,"/",new,"/g"," ", genome))
}

