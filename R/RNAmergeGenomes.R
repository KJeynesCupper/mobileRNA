#' Merge two FASTA genome assemblies
#'
#' @description Merges two reference genomes (.fa/.fasta). into one single
#' reference with modified chromosome names to ensure distinguishability.
#'
#'
#'@param genomeA path; directory path to a genome reference assembly file in
#'FASTA format (.fa/.fasta).
#'
#'@param genomeB path; directory path to a genome reference assembly file in
#'FASTA format (.fa/.fasta).
#'
#'@param out_dir path; a character string or a \code{base::connections()} open
#'for writing. Including file output name, and file extension of `.fa`. 
#'
#'@param abbreviationGenomeA character; string to represent prefix added to 
#'existing chromosome names in `genomeA`. Default set as "A", which is 
#'separated from existing chromosome names by an underscore (_). 
#'
#'@param abbreviationGenomeB character; string to represent prefix added to 
#'existing chromosome names in `genomeB`. Default set as "B", which is 
#'separated from existing chromosome names by an underscore (_). 
#'
#'@param cores logical;  the number of CPU cores for parallel computation. 
#'By default, `cores = TRUE` which tells the system to detect the number of CPU 
#'cores on the current host.While, `cores = FALSE` allows a user defined number
#'of cores which is stored in the variable `number_cores`. 
#'
#'@param number_cores numeric; number of CPU cores to use during parallel 
#'computation. 
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
#' # Load BiocFileCache package
#' library(BiocFileCache)
#' 
#' # generate cache in mobileRNA package
#' package_directory <- system.file(package = "mobileRNA")
#' 
#' # generate custom folder in package directory
#' cache_directory <- file.path(package_directory, "cache")
#' 
#' # Initialize a BiocFileCache in the cache directory
#' cache <-  BiocFileCache(cache_directory, ask = FALSE)
#' 
#' # Generate URL to FASTA files
#' url_remote <- "https://github.com/KJeynesCupper/assemblies/raw/main/"
#' fasta_1 <- paste0(url_remote, "chr12_Eggplant_V4.1.fa.gz")
#' fasta_2 <- paste0(url_remote,"chr2_S_lycopersicum_chromosomes.4.00.fa.gz")
#'  
#' # Add FASTA files to cache:
#' add_fasta1 <- bfcadd(cache,"chr12_Eggplant_V4.1", fpath=fasta_1)
#' rid3 <- names(add_fasta1)
#' 
#' add_fasta2 <- bfcadd(cache,"chr2_S_lycopersicum_chromosomes",fpath=fasta_2)
#' rid4 <- names(add_fasta2)
#' 
#' # rid3 and rid4 object contain the path location to the FASTA files 
#'
#' # run function to merge
#' merged_ref <- RNAmergeGenomes(genomeA = cache[[rid3]], genomeB = cache[[rid4]],
#'                                out_dir = tempfile("merged_genome",
#'                                tmpdir = cache_directory, 
#'                                fileext = ".fa"), 
#'                                cores = FALSE,number_cores = 1)
#'
#'
#'
#'
#' # or, to set specific changes to chromosome names:
#'
#' merged_ref2 <- RNAmergeGenomes(genomeA = cache[[rid3]], 
#'                                genomeB = cache[[rid4]] ,
#'             out_dir = tempfile("merged_genome2", tmpdir = cache_directory, 
#'             fileext = ".fa"),
#'             abbreviationGenomeA = "SM",
#'             abbreviationGenomeB = "SL",  
#'            cores = FALSE, number_cores = 1)
#'             
#' @importFrom Biostrings "readDNAStringSet"
#' @importFrom Biostrings "writeXStringSet"
#' @importFrom parallel "detectCores"
#' @importFrom parallel "makeCluster"
#' @importFrom parallel "parLapply"
#' @importFrom parallel "stopCluster"
#' @importFrom doParallel "registerDoParallel"
#' @importFrom foreach "foreach"
#' @importFrom foreach "%dopar%"
#' @importFrom parallel "makePSOCKcluster"
#' @importFrom BiocParallel MulticoreParam
#' @importFrom BiocParallel SnowParam
#' @importFrom BiocParallel bplapply
#' @export
RNAmergeGenomes <- function(genomeA, genomeB, out_dir,
                             abbreviationGenomeA = "A",
                             abbreviationGenomeB = "B",
                             cores = TRUE, number_cores = 2) {
  if (missing(genomeA)) {
    stop("Please specify genomeA, a connection to a FASTA file in local")
  }
  if (missing(genomeB)) {
    stop("Please specify annotationA, a connection to a FASTA file in local")
  }
  if (missing(out_dir) || !grepl("\\.fa$", out_dir)) {
    stop("Please specify out_dir, a connection to a local directory to write and 
          save merged annotation. Ensure file name with extension 
          (.fa or .fasta) is supplied.")
  }
  
  if (cores) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = number_cores)
  } else {
    BPPARAM <- BiocParallel::SnowParam(workers = number_cores)
  }
  
  ref1_fragments <- BiocParallel::bplapply(genomeA, function(file) {
    Biostrings::readDNAStringSet(file, format = "fasta")
  }, BPPARAM = BPPARAM)
  
  ref2_fragments <- BiocParallel::bplapply(genomeB, function(file) {
    Biostrings::readDNAStringSet(file, format = "fasta")
  }, BPPARAM = BPPARAM)
  
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
  message("Replacing chromosome names...")
  ref1_names <- names(ref1)
  ref1_newnames <- paste0(abbreviationGenomeA, "_", ref1_names)
  ref1_newnames <- sub("\\.", "", ref1_newnames)
  names(ref1) <- ref1_newnames
  
  ref2_names <- names(ref2)
  ref2_newnames <- paste0(abbreviationGenomeB, "_", ref2_names)
  ref2_newnames <- sub("\\.", "", ref2_newnames)
  names(ref2) <- ref2_newnames
  
  # Merge genomes
  merged_genome <- append(ref1, ref2)
  message("Genomes have been successfully merged.")
  message("Attempting to save merged genome to: ", out_dir, "\n")
  
  message("-- Note:Please be patient, this next step may take a long time \n")
  
  BiocParallel::bplapply(1:number_cores, function(i) {
    Biostrings::writeXStringSet(merged_genome,
                                file = paste0(dirname(out_dir), "/tempfile_", i,
                                              ".fa"),
                                format = "fasta")
  }, BPPARAM = BPPARAM)
  message("All temporary files have been created.")
  message("Merging temporary files... ")
  loc <- paste0(dirname(out_dir), "/tempfile_*.fa")
  system(paste0("cat ", loc, " > ", out_dir), intern = TRUE)
  
  message("Deleting temporary files...")
  system(paste0("rm ", loc), intern = TRUE)
  
  return(merged_genome)
}
