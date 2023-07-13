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
#'@param out_dir either a character string or a \code{base::connections()} open
#'for writing. Place path to output directory in "", including file output name.
#'Output name must have a ".fa" extension.
#'
#'@param abbreviationGenomeA a string placed in "", to replace chromosome names
#'within \code{genomeA}.Default set as "A".
#'
#'@param abbreviationGenomeB a string placed in "", to replace chromosome names
#'within \code{genomeB}.Default set as "B".
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
#'Please be aware that this function will take a very long time to process if 
#'you are working with large genomes which together take up more than 1Gb 
#'storage. Please be patient and allow the function to run. 
#' 
#' The functions primary goal is to merge two FASTA files, however, when
#' merging genomic files it is critical that the two genomes are distinguishable
#' by the chromosome names. As a default setting, the function extracts the
#' chromosome names for a given FASTA file and alters the name while retaining
#' the identifying number.
#'
#' The function requires the input of two reference genomes, where one
#' represents `Genome-A` and the other represents `Genome-B`. As default, the
#' function will rename the chromosome names in `Genome-A` to "A". For example,
#' A0, A1, A2 etc. To set a custom chromosome name for `Genome-A` alter the
#' argument \code{abbreviationGenomeA}. While, for  `Genome-B` as default the
#' chromosome names will be named "B", for example, B0, B1, B2 etc. To set a
#' custom chromosome name for `Genome-A` alter the argument
#' \code{abbreviationGenomeB}. The function can do so by  draw the chromosome
#' number within the given GFF file, remove all prior character or numerical
#' values, and replace it with the supplied string.
#'
#'Please note that this function uses parallel computation using the packages 
#'(`parallel`, `foreach` & `doParallel`) to improve the processing speed. Hence
#'common errors may occur is you have some parallel computing going on in the 
#'background. To improve the handling, the function will create a number of 
#'temporary files in your output directory, which are then merged into a single
#'file and deleted. 
#'
#'IMPORTANT:  The genome reference and annotation of a species
#'must have chromosomes with matching names. It is critical that if you used
#'the [mobileRNA::RNAmergeAnnotations()] function to to create a merged genome
#'annotation,that you treat the input references in the same way.
#'
#'
#'
#' @examples
#'
#'# URL for example genome assemblies
#'  url_remote <- "https://github.com/KJeynesCupper/assemblies/raw/main/"
#'
#'  fasta_1 <- paste0(url_remote, "chr12_Eggplant_V4.1.fa.gz")
#'
#'  fasta_2 <- paste0(url_remote,"chr2_S_lycopersicum_chromosomes.4.00.fa.gz")
#'
#'
#' # run function to merge
#' merged_ref <- RNAmergeGenomes(genomeA = fasta_1, genomeB = fasta_2,
#'                                out_dir = "../../merged_ref.fa", 
#'                                cores = FALSE,number_cores = 1)
#'
#'
#'
#'
#' ## or, to set specific changes to chromosome names:
#'  # genomeA represents Solanum melongena and the chromosomes will be
#'  # abbreviated to `SM.
#' # genomeB represents Solanum lycopersicum and the chromosomes will be
#' # abbreviated to `SL`.
#'
#' merged_ref2 <- RNAmergeGenomes(genomeA = fasta_1,
#'             genomeB = fasta_2 ,
#'             out_dir = "./merged_ref.fa",
#'             abbreviationGenomeA = "SM",
#'             abbreviationGenomeB = "SL",  
#'            cores = FALSE, number_cores = 1)
#'             
#'             
#'              
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
#' @export

RNAmergeGenomes <-  function(genomeA, genomeB, 
                             out_dir,abbreviationGenomeA = "A",
                             abbreviationGenomeB = "B", cores = TRUE, 
                             number_cores = 4) {
    
    if (base::missing(genomeA)) {
      stop("Please specify genomeA, a connection to a FASTA file in local")
    }
    
    if (base::missing(genomeB)) {
      stop("Please specify annotationA, a connection to a FASTA file in local")
    }
    
    if (base::missing(out_dir) || !grepl("\\.fa$", out_dir)) {
      stop("Please specify out_dir, a connection to a local directory to write 
           and save merged annotation. Ensure file name with extension 
           (.fa or .fasta) is supplied.")
    }
    # Load FASTA files in parallel
  if(cores){
    number_cores <- parallel::detectCores()
    cluster <- parallel::makeCluster(number_cores)
  }
  else if(cores == FALSE){
    number_cores <- number_cores
    cluster <- parallel::makePSOCKcluster(number_cores)
  }
    suppressWarnings(
      ref1_fragments <- parallel::parLapply(cluster, genomeA, function(file) {
        Biostrings::readDNAStringSet(file, format = "fasta")}))
    
    suppressWarnings(
      ref2_fragments <- parallel::parLapply(cluster, genomeB, function(file) {
        Biostrings::readDNAStringSet(file, format = "fasta") }))
    # Stop the cluster
    parallel::stopCluster(cluster)
    
    # join fragments together 
    ref1 <- do.call(c, ref1_fragments)
    ref2 <-do.call(c, ref2_fragments)
    
    # Replace chromosome names in reference genomes
    message("Replacing chromosome names")
    ref1_names <- names(ref1)
    ref1_newnames <- paste0(abbreviationGenomeA, "_", ref1_names)
    ref1_newnames <- sub("\\.", "", ref1_newnames)
    names(ref1) <- ref1_newnames
    
    ref2_names <- names(ref2)
    ref2_newnames <- paste0(abbreviationGenomeB, "_", ref2_names)
    ref2_newnames <- sub("\\.", "", ref2_newnames)
    names(ref2) <- ref2_newnames
    
    # merge genomes 
    merged_genome <- append(ref1, ref2)
    message("Genomes have been successfully merged")
    message("Attempting to save merged genome to: ", out_dir)
    
    `%dopar%` <- foreach::`%dopar%`
    # set number of cores for parelle, TRUE sets most for systm.  
    if(cores){
      number_cores <- parallel::detectCores()
      cluster <- parallel::makeCluster(number_cores)
    }
    else if(cores == FALSE){
      number_cores <- number_cores
      cluster <- parallel::makePSOCKcluster(number_cores)
      
    }
    
    doParallel::registerDoParallel(cluster)
    
    message("Please be patient, this next step may take a long time .. ")
    # Create a progress bar

    foreach::foreach(i = 1:number_cores) %dopar% {
      Biostrings::writeXStringSet(merged_genome,
                                  file = 
                                    paste0(dirname(out_dir), 
                                           "/tempfile_", i, ".fasta"), 
                                  format = "fasta")
    
    }
  
    on.exit(parallel::stopCluster(cluster))
    message("All temporary files have been created")
    message("Merging temporary files ... ")
    message("At this point get a cup of tea, its going to be a while ... ")
    
    loc <- paste0(dirname(out_dir), "/tempfile_*")
    system(paste0("cat ", loc, " > ", out_dir), intern = TRUE)
    
    message("Deleting temporary files ...")
    system(paste0("rm ", loc), intern = TRUE)
    return(merged_genome)
  }
  
