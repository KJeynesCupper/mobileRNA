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
#' As default, the function as a prefix to chromosomes which contains a either
#' "A" or "B" and "_". It will rename the chromosome names in `annotationA`
#' to "A_". For example, A_0, A_1, A_2 etc. To set a custom chromosome name for
#' `annotationA` alter the argument \code{abbreviationAnnoA}. While, for
#' `annotationB` as default the chromosome names will be named "B_", 
#' for example, B_0, B_1, B_2 etc. To set a custom chromosome name for 
#' `annotationB` alter the argument \code{abbreviationAnnoB}. Please note that 
#' the underscore is added automatically, hence, when setting a custom prefix 
#' just include character values. 
#' 
#' 
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
#'*IMPORTANT:* The genome reference and annotation of a species
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
#'@param annotationA path; directory path to genome annotation assembly file in
#'GFF format.
#'
#'
#'@param annotationB path; directory path to genome annotation assembly file in
#'GFF format.
#'
#'@param out_dir path; a character string or a \code{base::connections()} open
#'for writing. Including file output name, and file extension of `.fa`. 
#'
#'
#'@param abbreviationAnnoA character; string to represent prefix added to 
#'existing chromosome names in `annotationA`. Default set as "A", which is 
#'separated from existing chromosome names by an underscore (_). 
#'
#'
#'
#'@param abbreviationAnnoB character; string to represent prefix added to 
#'existing chromosome names in `annotationB`. Default set as "B", which is 
#'separated from existing chromosome names by an underscore (_). 
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
#' anno1 <- paste0(url_remote,"chr12_Eggplant_V4.1_function_IPR_final.gff")
#'
#' anno2 <- paste0(url_remote, "chr2_ITAG4.0_gene_models.gff")
#'
#'# edit and merge
#' merged_anno <- RNAmergeAnnotations(annotationA = anno1,annotationB = anno2,
#'              out_dir = tempfile("merged_annotation",tempdir(), 
#'              fileext = ".gff3"))
#'
#'
#'
#'
#'
#' @importFrom GenomeInfoDb "seqlevels"
#' @importFrom rtracklayer "export"
#' @importFrom GenomicRanges "GRangesList"
#' @importFrom rtracklayer "import"
#' @export
#'
RNAmergeAnnotations <- function(annotationA, annotationB,
                               out_dir,
                               abbreviationAnnoA = "A",
                               abbreviationAnnoB = "B", 
                               exportFormat = "GFF"){

  if (base::missing(annotationA) || !base::inherits(annotationA, "character")) {
    stop("Please specify annotationA object; path to GFF file.")
  }
  if (base::missing(annotationB) || !base::inherits(annotationB, "character")){
    stop("Please specify annotationA; path to GFF file.")
  }
  if (base::missing(out_dir) || grepl("\\.gff$", out_dir)) {
    stop("Please specify out_dir, a connection to a local directory to write and 
          save merged annotation. \n Ensure file name with extension (GFF) is 
          supplied.")
  }
  annotationA <- rtracklayer::import(annotationA)
  annotationB <- rtracklayer::import(annotationB)  
  message("Adding abbreviations to chomosome names ...  \n")
  # replace names with prefix and remove punc.
  annotationA_seqnames <- gsub("\\.", "", paste0(abbreviationAnnoA, "_",   
                                          GenomeInfoDb::seqlevels(annotationA)))
  annotationB_seqnames <- gsub("\\.", "", paste0(abbreviationAnnoB, "_",
                                          GenomeInfoDb::seqlevels(annotationB)))
  # rename 
  GenomeInfoDb::seqlevels(annotationA) <- annotationA_seqnames
  GenomeInfoDb::seqlevels(annotationB) <- annotationB_seqnames
  # Write out altered GFF files
  location <- dirname(out_dir)
  annoA_save <- paste0(location,"/", "annotationA_altered.gff")
  annoB_save <- paste0(location,"/", "annotationB_altered.gff")
  rtracklayer::export(annotationA, annoA_save, format = exportFormat)
  message("New annotation file annotationB: ", annoA_save, "\n")
  rtracklayer::export(annotationB, annoB_save, format = exportFormat)
  message("New annotation file created: ", annoB_save, "\n")
  message("Merging altered annotation files ... \n")
  gr_list <- GenomicRanges::GRangesList(annotationA, annotationB)
  concatenated_gff <- unlist(gr_list)
  rtracklayer::export(concatenated_gff, out_dir, format = exportFormat)
  message("New merged annnotation with modified chromosome names has been 
created saved to: ", out_dir)
  return(concatenated_gff)
}
