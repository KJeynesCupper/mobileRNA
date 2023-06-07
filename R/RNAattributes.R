#' Identify specific attributes associated to a sRNA cluster
#'
#' Based on genomic coordinates, assign sRNA clusters with an annotation that
#' has an exact match based on chromosome number, start and end coordinates
#'
#' @details The function merges an annotation (.gff/.gff3) file with
#' the sRNA data set based on the chromosome, start and end coordinates. It is
#' important that any alteration which were made to the genome reference
#' (FASTA), such as alterations to the chromosome name, must be carried forth
#' to the genome annotation file. If alterations were made to the reference
#' genome using the [mobileRNA::RNAmergeGenomes()] function, alteration
#' inline with these can be accomplished using the
#' [mobileRNA::RNAmergeAnnotations()] function.
#'
#' @param data data frame; containing rows of potential dicer-derived clusters
#' including columns which supply the genomic coordinates, where `chr` supplies
#' the chromosome number, `start` and `end` which supply the coordinates.
#'
#' @param annotation A path, URL, connection or GFFFile object. A genome
#' reference annotation file (.gff/.gff1/.gff2/.gff3).
#'
#'
#'@return Adds an additional 6 columns containing information which overlaps
#'any sRNA cluster loci. These columns represent the standard columns
#'in a GFF file.
#' @export
#' \dontrun{
#'
#' # find all overlapping attributes using a merged reference to the starting
#' # data
#' data("sRNA_data")
#'
#' attributes_df <- RNAattributes(data = sRNA_data,
#' annotation = "./annotation/merged/merged_annotation.gff3")
#'
#'
#' # find overlaps & annotate potential mobile sRNAs with overlapping regions
#' # in genome of destination tissues. The simulated data aims to simulate RNA
#' # movement from eggplant to tomato.
#' data("sRNA_data_mobile")
#' attributes_mobile_destination <- RNAattributes(data = sRNA_data_mobile,
#' annotation = "./annotation/tomato_annotation.gff")
#'
#' # find overlaps & annotate potential mobile sRNAs with overlapping regions in
#' # genome of origin to provide insight.
#' attributes_mobile_origin <- RNAattributes(data = sRNA_data_mobile,
#' annotation = "./annotation/eggplant_annotation.gff")
#'
#' }
#'
RNAattributes <- function(data, annotation){
  anno_data <- .import_annotations(annotation)
  res <- merge(data,anno_data, by=c("chr","start", "end"),all.x=TRUE)
  return(res)
}

