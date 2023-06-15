#' Merge two GFF genome annotation files
#'
#' @description Merges two genomes annotation files (`GFF format`) into one
#' single GFF file saved to desired directory. Adds prefix to each chromosome
#' name in the genome annotations to ensure each assembly is distinguishable
#' within the merged file, and saved the individually altered references to
#' local directory.
#'
#' @details
#' The functions primary goal is to merge two GFF files, however, when
#' merging genomic files it is critical that the two genomes are distinguishable
#' by the chromosome names. As a default setting, the function extracts the
#' chromosome names for a given GFF file and adds an ID before each chromosome
#' name and removes any periods.
#'
#'
#' As default, the function will rename the chromosome names in `annotationA`
#' to "A". For example, A0, A1, A2 etc. To set a custom chromosome name for
#' `annotationA` alter the argument \code{abbreviationAnnoA}. While, for
#' `annotationB` as default the chromosome names will be named "B", for example,
#' B0, B1, B2 etc. To set a custom chromosome name for `annotationB` alter the
#' argument \code{abbreviationAnnoB}.
#'
#' The merged genome is saved to the specified output directory, and requires
#' the user to set the name with a `gff` extension. The user must load the
#' annotations into the R Global Environment, hence, the alterations made by
#' this function will alter the variable. The altered annotations will also
#' be saved to the same directory as supplied for the output file, with a fixed
#' output name of `annotationA_altered.gff` and `annotationB_altered.gff`,
#' respectively.
#'
#'
#'
#'IMPORTANT:  The genome reference and annotation of a species
#'must have chromosomes with matching names. It is critical that if you used
#'the [mobileRNA::RNAmergeGenomes()] function to to create a merged reference
#'genome,that you treat the input annotations in the same way.
#'
#'@return
#'The function output the individually altered annotation files, plus, a merged
#'annotation file to the global environment and the given directory.
#'The altered annotation files contain either default or custom alterations to
#'the chromosome names and are written to the same directory as the merged file,
#' with the given names addition of `annotationA_altered.gff` and
#' `annotationB_altered.gff`, respectively.
#'
#'
#'
#'@param annotationA GRanges object; a genome annotation assembly file in
#'GFF format.
#'
#'
#'@param annotationB  GRanges object; a genome annotation assembly file in
#'GFF format.
#'
#'@param out_dir either a character string or a \code{base::connections()} open
#'for writing. Place path to output directory in "", including file output name
#'with extension of either `.gff` or `.gff3`.
#'
#'
#'@param abbreviationAnnoA a string placed in "", to replace chromosome names
#'within \code{annotationA}.Default set as "A".
#'
#'
#'
#'@param abbreviationAnnoB a string placed in "", to replace chromosome names
#'within \code{annotationB}.Default set as "B".
#'
#'
#'
#' @examples
#'
#' # import GFF files into R
#'
#'# Read reference genomes
#' url_remote <- "https://github.com/KJeynesCupper/assemblies/raw/main/"
#'
#' import_anno1_url <- paste0(url_remote, "chr12_Eggplant_V4.1_function_IPR_final.gff")
#'
#' import_anno2_url <- paste0(url_remote, "chr2_ITAG4.0_gene_models.gff")
#'
#'  # Read in GFF files
#' anno1 <- rtracklayer::import(import_anno1_url)
#' anno2 <- rtracklayer::import(import_anno2_url)
#'
#'# edit and merge
#' merged_anno <- RNAmergeAnnotations(annotationA = anno1,
#'                                  annotationB = anno2,
#'              out_dir = "../references/merged_annotation.gff3")
#'
#' ## Set specific pre-fixes:
#' ###  annotationA represents Solanum melongena, so abbreviated to `SM`
#' ###  annotationB represents Solanum lycopersicum, so abbreviated to `SL`
#'
#' merged_anno_2 <- RNAmergeAnnotations(annotationA = anno1,annotationB = anno2,
#'               out_dir = "../references/merged_annotation.gff3",
#'               abbreviationAnnoA = "SM",abbreviationAnnoB = "SL")
#'
#'
#'
#'
#'
#' @importFrom Repitools "annoGR2DF"
#' @importFrom utils "write.table"
#' @importFrom dplyr "bind_rows"
#' @export
#'



RNAmergeAnnotations <- function(annotationA, annotationB,
                               out_dir,
                               abbreviationAnnoA = "A",
                               abbreviationAnnoB = "B"){

  if (base::missing(annotationA) || !base::inherits(annotationA, "GRanges")) {
    stop(paste("Please specify annotationA object; GRanges object "))
  }

  if (base::missing(annotationB) || !base::inherits(annotationB, "GRanges")){
    stop(paste("Please specify annotationA; GRanges object"))
  }

  if (base::missing(out_dir) || grepl("\\.gff$", out_dir)) {
    stop(paste("Please specify out_dir, a connection to a local directory to
               write and save merged annotation. Ensure file name with extension
               (.gff or .gff3) is supplied."))
  }
  message("Adding abbreviations to chomosome names ... ")

  # Convert granges to dataframe
    annotationA <- Repitools::annoGR2DF(annotationA)
    annotationB <- Repitools::annoGR2DF(annotationB)
  # replace names with prefix and remove punc.
    annotationA$chr <- paste0(abbreviationAnnoA,"_", annotationA$chr)
    annotationA$chr <- sub("\\.", "", annotationA$chr)

    annotationB$chr <- paste0(abbreviationAnnoB,"_", annotationB$chr)
    annotationB$chr <- sub("\\.", "", annotationB$chr)

    # Write out altered GFF files
    location <- dirname(out_dir)
    annoA_save <- paste0(location,"/", "annotationA_altered.gff")
    annoB_save <- paste0(location,"/", "annotationB_altered.gff")

      utils::write.table(annotationA, file = annoA_save,sep = "\t",
                         quote = FALSE,row.names = FALSE, col.names = FALSE)

      message("New annotation file created: ", annoA_save)

      utils::write.table(annotationB, file = annoB_save,sep = "\t",
                         quote = FALSE, row.names = FALSE, col.names = FALSE)


    message("New annotation file created: ", annoB_save)
    message("Merging altered annotation files ...")

    # merge gff files
    concatenated_gff <- dplyr::bind_rows(annotationA, annotationB)
    utils::write.table(concatenated_gff, file = out_dir,sep = "\t",
                       quote = FALSE, row.names = FALSE, col.names = FALSE)

    return(concatenated_gff)

    message("New merged annnotation with modified chromosome names has been
    created saved to:", out_dir)
}
