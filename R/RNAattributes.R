#' Identify specific attributes associated to a sRNA cluster
#'
#' Based on genomic coordinates, assign sRNA clusters with an annotation.
#'
#' @details The function merges the annotation (.gff/.gff3) information with
#' the sRNA data set based on the chromosome, start and end coordinates.
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

