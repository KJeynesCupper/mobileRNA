#' Identify specific attributes associated to a sRNA cluster
#'
#' Based on genomic coordinates, assign sRNA clusters with an annotation that
#' has an exact match based on chromosome number, start and end coordinates
#'
#' @details The function merges an annotation (.gff/.gff3) file with
#' the sRNA data set based on the chromosome, start and end coordinates. It is
#' important that any alteration which were made to the genome reference (FASTA), such
#' as alterations to the chromosome name, must be carried forth to the
#' genome annotation file. If alterations were made to the reference genome
#' using the [RNAlocate::RNAmergeGenomes()] function, alteration inline with
#' these can be accomplished using the [RNAlocate::RNAmergeAnnotations()]
#' function.
#'
#' @param data data frame; containing rows of potential dicer-derived clusters
#' including columns which supply the genomic coordinates, where `chr` supplies the
#' chromosome number, `start` and `end` which supply the coordinates.
#'
#' @param annotation A path, URL, connection or GFFFile object. A genome
#' reference annotation file (.gff/.gff1/.gff2/.gff3).
#'
RNAattributes <- function(data, annotation){
  anno_data <- .import_annotations(annotation)
  res <- merge(data,anno_data, by=c("chr","start", "end"),all.x=TRUE)
  return(res)
}

