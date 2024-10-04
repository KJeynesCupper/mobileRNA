#' Merge two genome reference assemblies (FASTA format)
#'
#' @description Merges two reference genomes (FASTA) into one single
#' reference with modified chromosome names. 
#' 
#' Typically, use genomeA as the origin tissue genome assembly, and genomeB as 
#' the genome from which mobile RNAs are produced by. 
#'
#'@param genomeA path; path to a genome reference assembly file in FASTA format.
#'
#'@param genomeB path; path to a genome reference assembly file in FASTA format.
#'
#'@param output_file path; a character string or a \code{base::connections()} 
#'open for writing. Including file output name, and file extension of `.fa` or 
#'`.fasta`. 
#'
#'@param GenomeA.ID character; string to represent prefix added to 
#'existing chromosome names in `genomeA`. Default set as "A"
#'
#'@param GenomeB.ID character; string to represent prefix added to 
#'existing chromosome names in `genomeB`. Default set as "B". 
#'
#'@param compress.output logical; state whether the output file should be in 
#'a compressed format (gzip)
#'
#'@return Returns a single FASTA format file containing both genome assemblies
#'with edited chromosome names (prefixes, and removal of periods) to the given
#'directory. 
#'
#'@details
#' The function merges two FASTA files, however, when merging genomic files it 
#' is critical that the two genomes are distinguishable by the chromosome names. 
#' As a default setting, the function extracts the chromosome names for the
#' given FASTA files and alters adds a unique prefix while retaining the 
#' identifying number. Plus, removes any periods. 
#'
#' As default, the function will rename the chromosome names in `genome_A` to 
#' "A" and separates the prefix and the existing chromosome names with an 
#' underscore ("_"). For example, A_0, A_1, A_2 etc.
#' 
#' Please note that the underscore is added automatically, hence, when setting a
#' custom prefix just includes character values. 
#' 
#'
#'**IMPORTANT:**  The genome reference and annotation of the same 
#'species/accession/variety must have chromosomes with matching names. It is 
#'critical that if you use the [mobileRNA::RNAmergeAnnotations()] function to
#' create a merged genome annotation, that you treat the input references in the 
#' same way.
#'
#'
#' @examples
#' fasta_1 <- system.file("extdata","reduced_chr12_Eggplant.fa.gz", 
#' package="mobileRNA")
#' 
#' fasta_2 <-system.file("extdata","reduced_chr2_Tomato.fa.gz",
#' package="mobileRNA")
#' 
#' 
#' output_file <- file.path(tempfile("merged_annotation", fileext = ".fa"))
#' 
#' merged_ref <- RNAmergeGenomes(genomeA = fasta_1, 
#' genomeB = fasta_2,
#' output_file = output_file)
#' 
#' 
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings writeXStringSet
#' @importFrom progress progress_bar
#' @export
RNAmergeGenomes <- function(genomeA, genomeB, output_file,
                             GenomeA.ID = "A",
                             GenomeB.ID = "B", 
                            compress.output = FALSE) {
  if (missing(genomeA) || !file.exists(genomeA)) {
    stop("Please specify genomeA, a connection to a FASTA file in local")
  }
  if (missing(genomeB) || !file.exists(genomeB)) {
    stop("Please specify annotationA, a connection to a FASTA file in local")
  }
  if (missing(output_file) || !grepl("\\.(fa|fasta|fasta.gz|fa.gz|fsa|fsa.gz)$", 
                                     output_file)) {
    stop("Please specify output_file, a connection to a local directory to write 
    and save merged annotation. Ensure file name with extension (.fa or .fasta) 
         is supplied.")
  }
  
  if (compress.output == TRUE || grepl("\\.gz$", output_file)) {
    stop("If you wish to compress the output file, please ensure the file 
      extension includes '.gz' ")
  }
  message("-- Note: Please be patient, this next step may take a long time --")

  # progress bar
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent :current/:total :eta",
    total = 5)
  pb$tick(0)  # start the progress bar
    # 1- load each genome into R 
  ref1 <-  Biostrings::readDNAStringSet(genomeA)

  pb$tick()
  ref2 <-  Biostrings::readDNAStringSet(genomeB)
    
  pb$tick()
  
  
  # Replace chromosome names in reference genomes
  ref1_names <- names(ref1)
  ref1_newnames <- paste0(GenomeA.ID, "_", ref1_names)
  ref1_newnames <- sub("\\.", "", ref1_newnames)
  names(ref1) <- ref1_newnames
  
  ref2_names <- names(ref2)
  ref2_newnames <- paste0(GenomeB.ID, "_", ref2_names)
  ref2_newnames <- sub("\\.", "", ref2_newnames)
  names(ref2) <- ref2_newnames
  pb$tick()
  # Merge genomes
  merged_genome <- append(ref1, ref2)
  pb$tick()
  # save merged genome
  Biostrings::writeXStringSet(merged_genome, 
                              output_file, 
                              format="fasta",
                              compress = compress.output)
  pb$tick()
  message("Output file has been saved to: ", output_file)

  return(merged_genome)
}
