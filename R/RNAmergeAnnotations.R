#' Merge two GFF genome annotation files
#'
#' @description Merges two genomes annotation files (`GFF format`) into one
#' single GFF file saved to desired directory.
#'
#' @details
#' The functions primary goal is to merge two GFF files, however, when
#' merging genomic files it is critical that the two genomes are distinguishable
#' by the chromosome names. As a default setting, the function extracts the
#' chromosome names for a given GFF file and adds an ID before each chromosome
#' name and removes any periods.
#'
#'
#'
#' The function requires the input of two reference annotation genomes, where one represents
#' `annotationA` and the other represents `annotationB`. As default, the function
#' will rename the chromosome names in `annotationA` to "A". For example, A0, A1 ,
#' A2 etc. To set a custom chromosome name for `annotationA` alter the argument
#' \code{abbreviationAnnoA}. While, for  `annotationB` as default the chromosome
#' names will be named "B", for example, B0, B1, B2 etc. To set a custom
#' chromosome name for `annotationB` alter the argument \code{abbreviationAnnoB}.
#' The function can do so by  draw the chromosome number within the given GFF file,
#' remove all prior character or numerical values, and replace it with the
#' supplied string.
#'
#'IMPORTANT:  The genome reference and annotation of a species
#'must have chromosomes with matching names. It is critical that if you used
#'the [mobileRNA::RNAmergeGenomes()] function to to create a merged reference genome,
#'that you treat the input annotations in the same way.
#'
#'@return
#'The function output the individually altered annotation files, plus, a merged
#'annotation file made up of the two altered files. The altered annotation files
#'contain either default or custom alterations to the chromosome names and are
#' written to the original file connection with the addition of "_altered.gff3"
#' at the end of the file name.
#'
#'The  merged annotation file is outputted to the global object variable when
#'assigned, as well as writing the merged annotation file in
#'GFF format to a give output directory.
#'
#'
#'@param annotationA the path  \code{base::connection()} to a genome annotation
#'file in GFF3 format
#'
#'@param annotationB the path  \code{base::connection()} to a genome annotation
#'file in GFF3 format
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
#'@param replace_chr_names a logical value indicating whether the chromosome
#'names of the supplied annotation files are to be altered or not. As default,
#'chromosome names are altered.
#'
#'
#' @examples \dontrun{
#' merged_anno <- RNAannotationMerge(annotationA = "./workplace/reference/ref1_annotation.gff3",
#'             annotationB = "./workplace/reference/ref2_annotation.gff3",
#'             out_dir = "./workplace/reference/merge/merged_annotation.gff3")
#'
#' ## or, to set specific changes to chromosome names. annotationA represents
#' ## the Solanum lycopersicum and the chromosomes will be abbreviated to `SL`,
#' ## and annotationB represents Solanum melongena and the chromosomes will be
#' ## abbreviated to `SM`.
#'
#' merged_anno_2 <- RNAannotationMerge(annotationA = "./workplace/reference/ref1_annotation.gff3",
#'             annotationB = "./workplace/reference/ref2_annotation.gff3",
#'             out_dir = "./workplace/reference/merge/merged_annotation.gff3",
#'             abbreviationAnnoA = "SL",
#'             abbreviationAnnoB = "SM")
#'
#'
#'
#'}
#'
#'
#' @importFrom utils "read.table"
#' @importFrom utils "write.table"
#' @importFrom tools "file_ext"
#' @export
RNAmergeAnnotations <- function(annotationA, annotationB,
                               out_dir,
                               abbreviationAnnoA = "A",
                               abbreviationAnnoB = "B",
                               replace_chr_names = TRUE){

  if (base::missing(annotationA) || !base::inherits(annotationA, c("character"))) {
    stop(paste("Please specify annotationA, a connection to a GFF3 file in local"))
  }

  if (base::missing(annotationB) || !base::inherits(annotationB, c("character"))) {
    stop(paste("Please specify annotationA, a connection to a GFF3 file in local"))
  }

  if (base::missing(out_dir) || !base::inherits(out_dir, c("character")) || tools::file_ext(out_dir == c("gff", "gff3"))) {
    stop(paste("Please specify out_dir, a connection to a local directory to
               write and save merged annotation. Ensure file name with extension
               (.gff or .gff3) is supplied."))
  }


  # Read in GFF files
  genome1 <- utils::read.table(annotationA, sep = "\t", fill = TRUE, comment.char = "#")
  genome2 <- utils::read.table(annotationB, sep = "\t", fill = TRUE, comment.char = "#")

  if (replace_chr_names) {
    message("Adding abbreviations to chomosome names ... ")
    # Replace chromosome names in both
    genome1$V1 <- paste0(abbreviationAnnoA,"_", genome1$V1 )
    genome2$V1 <- paste0(abbreviationAnnoB,"_", genome2$V1 )

    # remove any dots/periods in names
    genome1$V1 <- sub("\\.", "", genome1$V1)
    genome2$V1 <- sub("\\.", "", genome2$V1)

    # Write out altered GFF files
      utils::write.table(genome1, file = paste0(gsub("\\.gff3*", "", annotationA), "_altered.gff3"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      message("New annotation file created:  " , paste0(gsub("\\.gff3", "", annotationA), "_altered.gff3"))

      utils::write.table(genome2, file = paste0(gsub("\\.gff3*", "", annotationB), "_altered.gff3"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      message("New annotation file created:  " , paste0(gsub("\\.gff3", "", annotationB), "_altered.gff3"))

    # select files, and merge and save to new file
    file_1 <- paste0(gsub("\\.gff3*", "", annotationA), "_altered.gff3")
    file_2 <- paste0(gsub("\\.gff3*", "", annotationB), "_altered.gff3")
    message("Merging altered annotation files ...")
    system(paste0("cat ", file_1, " ", file_2, " >", out_dir))

    } else
      if (replace_chr_names == FALSE) {
         message("Merging un-altered annotation files ...")
        # select files, and merge and save to new file
        file_1 <- paste0(gsub("\\.gff3*", "", annotationA), "_altered.gff3")
        file_2 <- paste0(gsub("\\.gff3*", "", annotationB), "_altered.gff3")

        system(paste0("cat ", file_1, " ", file_2, " >", out_dir))
            }
    # Return concatenated GFF data.frame
    concatenated_gff <- utils::read.table(out_dir, sep = "\t", fill = TRUE, comment.char = "#")
    return(concatenated_gff)

    message("New merged annnotation with modified chromosome names has been created
            and saved to:", out_dir)
  }
