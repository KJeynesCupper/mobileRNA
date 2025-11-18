#' Merge two genome annotation files (GFF Format)
#'
#' @description Merges two genomes annotation files (GFF) into one
#' single GFF format file saved to the desired directory. This function adds a
#' unique prefix to the chromosome names in each genome annotation to ensure
#' each is distinguishable within the merged file.
#'
#' @details
#' As default, the function removes periods and adds a prefix to the existing
#' chromosome names. The prefix is separated from the original chromosome name
#' by an underscore. For example, based on the default settings, it will add the
#' prefix "A_" to the chromosome names in `annotationA`, for instance,
#'  A_0, A_1, A_2 etc.
#'
#' The merged genome is saved to the specified output directory, and requires
#' the user to set the name with a GFF format.
#'
#'
#'
#'**IMPORTANT:** The genome reference and annotation of a species
#'must have chromosomes with matching names. It is critical that if you used
#'the [mobileRNA::RNAmergeGenomes()] function to create a merged reference
#'genome,that you treat the input annotations in the same way.
#'
#'@return
#'A GFF format file containing the annotations of two genomes distinguishable by
#'the appended prefixes.
#'
#'
#'@param annotationA path; path to a genome annotation assembly file in GFF format.
#'
#'@param annotationB path; path to a genome annotation assembly file in GFF format.
#'
#'@param output_file path; a character string or a \code{base::connections()}
#'open for writing. Including file output name, and must have a GFF file
#'extension.
#'
#'@param AnnoA.ID character; string to represent prefix added to
#'existing chromosome names in `annotationA`. Default set as "A".
#'
#'@param AnnoB.ID character; string to represent prefix added to
#'existing chromosome names in `annotationB`. Default set as "B".
#'
#'@param format format of GFF output, either "gff", "gff1", "gff2", "gff3."
#'  Default is "gff3".
#'
#'
#' @examples
#'
#' anno1 <- system.file("extdata", "reduced_chr12_Eggplant.gff.gz",
#' package="mobileRNA")
#'
#' anno2 <- system.file("extdata","reduced_chr2_Tomato.gff.gz",
#' package="mobileRNA")
#'
#' output_file <- tempfile("merged_annotation", fileext = ".gff3")
#'
#' merged_anno <- RNAmergeAnnotations(annotationA = anno1, annotationB = anno2,
#'                                    output_file = output_file)
#'
#' @importFrom rtracklayer import
#' @importFrom rtracklayer export
#' @importFrom GenomicRanges GRangesList
#' @importFrom progress progress_bar
#' @importFrom S4Vectors elementMetadata
#' @importFrom GenomeInfoDb seqlevels
#' @export
#'
RNAmergeAnnotations <- function(annotationA, annotationB,
                               output_file,
                               AnnoA.ID = "A",
                               AnnoB.ID = "B",
                               format = "gff3"){

  if (base::missing(annotationA) || !base::inherits(annotationA, "character") ||
      !file.exists(annotationA)) {
    stop("Please specify annotationA object; path to GFF file.")
  }
  if (base::missing(annotationB) || !base::inherits(annotationB, "character") ||
      !file.exists(annotationB)){
    stop("Please specify annotationA; path to GFF file.")
  }
  if (base::missing(output_file) ||
      !grepl("\\.(gff)$|\\.(gff1)$|\\.(gff2)$|\\.(gff3)$|\\.(gff.gz)$|
             \\.(gff1.gz)$|\\.(gff2.gz)$|\\.(gff3.gz)$",output_file)) {
    stop("Please specify output_file, a connection to a local directory to write
          and save merged annotation.
          Ensure file contains suitable GFF extension.")
  }
  # start the progress bar
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent :current/:total :eta",
    total = 4)
  pb$tick(0)

  # import gff
  annotationA <- rtracklayer::import(annotationA)
  annotationB <- rtracklayer::import(annotationB)
  pb$tick()

  # remove full-stop
  annotationA_seqnames <- gsub("\\.", "", paste0(AnnoA.ID, "_",
                                                 GenomeInfoDb::seqlevels(annotationA)))
  annotationB_seqnames <- gsub("\\.", "", paste0(AnnoB.ID, "_",
                                                 GenomeInfoDb::seqlevels(annotationB)))

  # change seqnames
  GenomeInfoDb::seqlevels(annotationA) <- annotationA_seqnames
  GenomeInfoDb::seqlevels(annotationB) <- annotationB_seqnames

  pb$tick()
 # change classes
  annotationA_gff <- convertChar2Factor(annotationA)
  annotationB_gff <- convertChar2Factor(annotationB)

  pb$tick()
  # join
  gr_list <- GenomicRanges::GRangesList(annotationA_gff, annotationB_gff)
  concatenated_gff <- unlist(gr_list)
  # export merged
  rtracklayer::export(concatenated_gff, output_file, format=format)
  pb$tick()
  message("\nOutput file has been saved to: ", output_file)

  return(concatenated_gff)
}
