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
#'@param output_file path; a character string or a \code{base::connections()} 
#'open for writing. Including file output name, and must have a `.gff3` 
#'extension. 
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
#' @param format format of GFF output, either gff, gff1, gff2, gff3. Default is 
#' gff3.  
#'
#'
#' @examples
#' 
#' anno1 <- system.file("extdata", "reduced_chr12_Eggplant.gff", 
#' package="mobileRNA")
#' 
#' anno2 <- system.file("extdata","reduced_chr2_Tomato.gff", package="mobileRNA")
#' 
#' # Merge annotations and write them to a file in output_file.
#' # For this example, the result is written to a temporary file.
#' # For real use cases, the output_file should be an appropriate location on 
#' # your computer.
#' 
#' output_file <- tempfile("merged_annotation", fileext = ".gff3")
#' 
#' merged_anno <- RNAmergeAnnotations(annotationA = anno1, annotationB = anno2,
#'                                    output_file = output_file)
#'                                    
#' @importFrom rtracklayer "import"
#' @importFrom rtracklayer "export"
#' @importFrom GenomicRanges "GRangesList"
#' @importFrom progress "progress_bar"
#' @importFrom S4Vectors "elementMetadata" 
#' @importFrom GenomeInfoDb "seqlevels"
#' @export
#'
RNAmergeAnnotations <- function(annotationA, annotationB,
                               output_file,
                               abbreviationAnnoA = "A",
                               abbreviationAnnoB = "B",
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
      !grepl("\\.(gff3)$",output_file)) {
    stop("Please specify output_file, a connection to a local directory to write 
    and save merged annotation. \n Ensure file name with .gff3 extension is 
    supplied.")
  }
  # start the progress bar
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent :current/:total :eta",
    total = 6)
  pb$tick(0)  
  
  # import gff
  annotationA <- rtracklayer::import(annotationA)
  annotationB <- rtracklayer::import(annotationB) 
  pb$tick()
  
  # remove full-stop
  annotationA_seqnames <- gsub("\\.", "", paste0(abbreviationAnnoA, "_",   
                                                 GenomeInfoDb::seqlevels(annotationA)))
  annotationB_seqnames <- gsub("\\.", "", paste0(abbreviationAnnoB, "_",
                                                 GenomeInfoDb::seqlevels(annotationB)))
  
  # change seqnames
  GenomeInfoDb::seqlevels(annotationA) <- annotationA_seqnames
  GenomeInfoDb::seqlevels(annotationB) <- annotationB_seqnames

  pb$tick()
 # change classes 
  annotationA_gff <- convertChar2Factor(annotationA)
  annotationB_gff <- convertChar2Factor(annotationB)
  # output location
  location <- dirname(output_file)
  annoA_save <- file.path(location, paste0("annotationA_altered.", format))
  annoB_save <- file.path(location, paste0("annotationB_altered.", format))
  pb$tick() 

  # export individual gffs
  rtracklayer::export(annotationA_gff, annoA_save)
  rtracklayer::export(annotationB_gff, annoB_save)
  pb$tick()
  # join 
  gr_list <- GenomicRanges::GRangesList(annotationA_gff, annotationB_gff)
  concatenated_gff <- unlist(gr_list)
  # export merged 
  rtracklayer::export(concatenated_gff, output_file, format=format)
  pb$tick()
  cat("\n")
  message("Output files have been saved to: ")
  message("---- annotationA:", annoA_save)
  message("---- annotationB:", annoB_save)
  message("---- Merged genome:", output_file)
 
  return(concatenated_gff)
}
