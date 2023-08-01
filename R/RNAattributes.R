#' Overlap genome annotation file information with sRNA/mRNA-seq data 
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
#' @param input character; must be either "sRNA" or "mRNA"
#'@return Adds an additional columns from GFF containing information which
#'overlaps any sRNA cluster loci. These columns represent the standard columns
#'in a GFF file.
#' @export
#' @importFrom rtracklayer "import"
#' @importFrom Repitools "annoGR2DF"
#' @importFrom S4Vectors "mcols"
#' @importFrom stats "complete.cases"
#' @examples
#' \dontrun{
#'
#' data("sRNA_data")
#'
#' attributes_df <- RNAattributes(data = sRNA_data,
#'                     annotation = "./annotation/merged_annotation.gff3")
#'
#'
#'
#' }
#'
RNAattributes <- function(data, annotation, input= c("sRNA", "mRNA")){
  if (base::missing(data)) {
    stop("data is missing. data must be an object of class matrix, data.frame, 
         DataFrame")
  }
  if (base::missing(annotation)) {
    stop("annotation is missing. annotation must be an object of GFF format.")
  }
  cat("Please be patient...")
  anno_data <- rtracklayer::import(annotation)
  conversion <- Repitools::annoGR2DF(anno_data)
  if(input == "sRNA"){
    # check chromosome names match:
    tryCatch(
      {  
        res <- merge(data,conversion, by=c("chr","start", "end"),all.x=TRUE)
      }, error = function(e) {
        cat("An error occurred .. :", conditionMessage(e), "\n")
        cat("An error occurred ... chromosome names.", "\n")
      }
    )
    } else
      if (input == "mRNA"){
        genes_info <- conversion[which(S4Vectors::mcols(conversion)$type == "gene")]
        merged_gene_info <- merge(data, genes_info, by = "Gene", all.x = TRUE)
        res <- merged_gene_info[stats::complete.cases(merged_gene_info), ]
      }
  return(res)
}

