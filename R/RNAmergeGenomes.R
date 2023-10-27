#' Merge two FASTA genome assemblies
#'
#' @description Merges two reference genomes (.fa/.fasta). into one single
#' reference with modified chromosome names to ensure distinguishability.
#' Typically, use genomeA as the origin tissue genome assembly, and genomeB as 
#' the mobile/foreign genome. 
#'
#'@param genomeA path; directory path to a genome reference assembly file in
#'FASTA format (.fa/.fasta).
#'
#'@param genomeB path; directory path to a genome reference assembly file in
#'FASTA format (.fa/.fasta).
#'
#'@param output_file path; a character string or a \code{base::connections()} 
#'open for writing. Including file output name, and file extension of `.fa`. 
#'
#'@param abbreviationGenomeA character; string to represent prefix added to 
#'existing chromosome names in `genomeA`. Default set as "A", which is 
#'separated from existing chromosome names by an underscore (_). 
#'
#'@param abbreviationGenomeB character; string to represent prefix added to 
#'existing chromosome names in `genomeB`. Default set as "B", which is 
#'separated from existing chromosome names by an underscore (_). 
#'
#'@param BPPARAM An optional [BiocParallelParam()] instance determining the 
#'parallel back-end to be used during evaluation, or a list of BiocParallelParam 
#'instances, to be applied in sequence for nested calls to *BiocParallel*
#'functions.
#'
#'
#'@return Returns a single FASTA format file containing both  genome assemblies
#'with edited chromosome names (prefixes, and removal of periods) to the give
#'directory. 
#'@details
#' The function merges two FASTA files, however, when merging genomic files it 
#' is critical that the two genomes are distinguishable by the chromosome names. 
#' As a default setting, the function extracts the chromosome names for the
#' given FASTA files and alters adds a unique prefix while retaining the 
#' identifying number.
#'
#' The function requires the input of two FASTA reference genomes, where one
#' represents `genome_A` and the other represents `genome_B`. As default, the
#' function will rename the chromosome names in `genome_A` to "A_". For example,
#' A_0, A_1, A_2 etc. To set a custom chromosome name for `genome_A` alter the
#' argument \code{abbreviationGenomeA}. While, for  `genome_B` as default the
#' chromosome names will be named "B_", for example, B_0, B_1, B_2 etc. To set a
#' custom chromosome name for `genome_B` alter the argument
#' \code{abbreviationGenomeB}.  Please note that the underscore is added 
#' automatically, hence, when setting a custom prefix just include character 
#' values. 
#' 
#'Please be aware that this function will take a very long time to process if 
#'you are working with large genomes which together take up more than 1Gb 
#'storage. Please be patient and allow the function to run. 
#' 
#'Please note that this function uses parallel computation to improve the 
#'processing speed.
#'
#'*IMPORTANT:*  The genome reference and annotation of the same 
#'species/accession/variety must have chromosomes with matching names. It is 
#'critical that if you use the [mobileRNA::RNAmergeAnnotations()] function to
#' create a merged genome annotation,that you treat the input references in the 
#' same way.
#'
#'
#'
#' @examples
#' 
#' # Initialize a cache directory, using the BiocFileCache package, to store the 
#' # downloads used in this example
#' library(BiocFileCache)
#' cache_dir <- tools::R_user_dir("mobileRNA", which = "cache")
#' cache <- BiocFileCache(cache_dir)
#' 
#' 
#' # Construct URL to example FASTA files
#' url_remote <- "https://github.com/KJeynesCupper/assemblies/raw/main/"
#' 
#' fasta_1_url <- file.path(url_remote, "chr12_Eggplant_V4.1.fa.gz")
#' fasta_2_url <- file.path(url_remote,"chr2_S_lycopersicum_chromosomes.4.00.fa.gz")
#' 
#' # Download example FASTA files and add them to cache
#' fasta_1 <- bfcrpath(cache, fasta_1_url)
#' fasta_2 <- bfcrpath(cache, fasta_2_url)
#' 
#' 
#' # Merge FASTA files and write them to a file in output_file.
#' # For this example, the result is written to a temporary file.
#' # For real use cases, the output_file should be an appropriate location on your 
#' # computer.
#' 
#' output_file <- tempfile("merged_annotation", fileext = ".fa")
#' 
#' merged_ref <- RNAmergeGenomes(genomeA = fasta_1, 
#' genomeB = fasta_2,
#' output_file = output_file)
#'             
#' @importFrom Biostrings "readDNAStringSet"
#' @importFrom Biostrings "writeXStringSet"
#' @importFrom BiocParallel bplapply
#' @importFrom progress progress_bar
#' @export
RNAmergeGenomes <- function(genomeA, genomeB, output_file,
                             abbreviationGenomeA = "A",
                             abbreviationGenomeB = "B",
                            BPPARAM = BiocParallel::SerialParam()) {
  if (missing(genomeA) || !file.exists(genomeA)) {
    stop("Please specify genomeA, a connection to a FASTA file in local")
  }
  if (missing(genomeB) || !file.exists(genomeB)) {
    stop("Please specify annotationA, a connection to a FASTA file in local")
  }
  if (missing(output_file) || !grepl("\\.fa$", output_file)) {
    stop("Please specify output_file, a connection to a local directory to write 
    and save merged annotation. Ensure file name with extension (.fa or .fasta) 
         is supplied.")
  }
  message("-- Note: Please be patient, this next step may take a long time -- \n")

  # progress bar
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent :current/:total :eta",
    total = 5)
  pb$tick(0)  # start the progress bar
    # 1- load each genome into R 
  #ref1_fragments <- Biostrings::readDNAStringSet(genomeA)
  ref1_fragments <- BiocParallel::bplapply(genomeA, 
                                           FUN = Biostrings::readDNAStringSet,
                                           BPPARAM=BPPARAM)
  pb$tick()
  ref2_fragments <-  BiocParallel::bplapply(genomeB, 
                                            FUN = Biostrings::readDNAStringSet,
                                            BPPARAM=BPPARAM)
    
  pb$tick()
  if(!is.null(names(ref1_fragments))){
    ref1_fragments <- unname(ref1_fragments)
  }
  if(!is.null(names(ref2_fragments))){
    ref2_fragments <- unname(ref2_fragments)
  }
  
  # Join fragments together
  ref1 <- do.call(c, ref1_fragments)
  ref2 <- do.call(c, ref2_fragments)
  
  # Replace chromosome names in reference genomes
  ref1_names <- names(ref1)
  ref1_newnames <- paste0(abbreviationGenomeA, "_", ref1_names)
  ref1_newnames <- sub("\\.", "", ref1_newnames)
  names(ref1) <- ref1_newnames
  
  ref2_names <- names(ref2)
  ref2_newnames <- paste0(abbreviationGenomeB, "_", ref2_names)
  ref2_newnames <- sub("\\.", "", ref2_newnames)
  names(ref2) <- ref2_newnames
  pb$tick()
  # Merge genomes
  merged_genome <- append(ref1, ref2)
  pb$tick()
  # save merged genome
  Biostrings::writeXStringSet(merged_genome, 
                              output_file, 
                              format="fasta")
  pb$tick()
  message("Merged genome has been save to: ", output_file, "\n")
  
  return(merged_genome)
}
