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
#'@param exportFormat The format of the annotation output. If missing or the 
#'format does not match the file type, an error will occur. Default set to "GFF"
#'format. 
#'
#'
#' @examples
#' # import GFF files into R
#'
#'# Read reference genomes
#' url_remote <- "https://github.com/KJeynesCupper/assemblies/raw/main/"
#'
#' anno1_url <- paste0(url_remote,"chr12_Eggplant_V4.1_function_IPR_final.gff")
#'
#' anno2_url <- paste0(url_remote, "chr2_ITAG4.0_gene_models.gff")
#'
#'  # Read in GFF files
#' anno1 <- rtracklayer::import(anno1_url)
#' anno2 <- rtracklayer::import(anno2_url)
#'
#'
#'# edit and merge
#' merged_anno <- RNAmergeAnnotations(annotationA = anno1,annotationB = anno2,
#'              out_dir = "./merged_annotation.gff3")
#'
#'
#'
#'
#'
#' @importFrom GenomeInfoDb "seqlevels"
#' @importFrom rtracklayer "export"
#' @importFrom GenomicRanges "GRangesList"
#' @export
#'
RNAmergeAnnotations <- function(annotationA, annotationB,
                               out_dir,
                               abbreviationAnnoA = "A",
                               abbreviationAnnoB = "B", 
                               exportFormat = "GFF"){

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
  cat("Adding abbreviations to chomosome names ...  \n")
  # replace names with prefix and remove punc.
  annotationA_seqnames <- gsub("\\.", "", paste0(abbreviationAnnoA, "_",   GenomeInfoDb::seqlevels(annotationA)))
  annotationB_seqnames <- gsub("\\.", "", paste0(abbreviationAnnoB, "_",   GenomeInfoDb::seqlevels(annotationB)))
  
  # rename 
  GenomeInfoDb::seqlevels(annotationA) <- annotationA_seqnames
  GenomeInfoDb::seqlevels(annotationB) <- annotationB_seqnames
  
  # Write out altered GFF files
  location <- dirname(out_dir)
  annoA_save <- paste0(location,"/", "annotationA_altered.gff")
  annoB_save <- paste0(location,"/", "annotationB_altered.gff")
  
  rtracklayer::export(annotationA, annoA_save, format = exportFormat)
  cat("New annotation file annotationB: ", annoA_save, "\n")
  rtracklayer::export(annotationB, annoB_save, format = exportFormat)
  cat("New annotation file created: ", annoB_save, "\n")
  
  cat("Merging altered annotation files ... \n")
  gr_list <- GenomicRanges::GRangesList(annotationA, annotationB)
  concatenated_gff <- unlist(gr_list)
  
  rtracklayer::export(concatenated_gff, out_dir, format = exportFormat)
  cat("New merged annnotation with modified chromosome names has been created saved to: ", out_dir)
  return(concatenated_gff)
}
