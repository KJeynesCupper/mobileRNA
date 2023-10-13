#' Map & cluster sRNAseq reads using `ShortStack`
#'
#' @description Map sRNAseq FASTQ samples to a given genome assembly using 
#' ShortStack. 
#' 
#' @details
#' This function undertakes our recommended workflow for mapping sRNAseq reads
#' to a given genome for the analysis of small RNAs in a biological system. 
#' See appendix of vignette for manual pipeline. 
#' 
#' The `sRNAmapper()` function invokes a number of OS commands, and is dependent 
#' on the installation of `ShortStack` (>= 4.0) with Conda.  Please refer to 
#' \link{https://github.com/MikeAxtell/ShortStack} for more information. Please 
#' note that `ShortStack` is only compatible with Linux and Mac operating 
#' systems. 
#' 
#' The pipeline undertakes de novo detection of sRNA-producing loci and
#' alignment, where the output of each are stored in their respective folders in 
#' the users desired location. The de novo detection of sRNA-producing loci
#' analyses each sample to identify de novo sRNA-producing loci (ie. sRNA
#' clusters), and joins these results into a single file called "locifile.txt". 
#' The alignment step aligns and clusters each sample to the genome reference
#' along with the file containing the de novo sRNA clusters. 
#' 
#' @references ShortStack \link{https://github.com/MikeAxtell/ShortStack}
#' @param input_files_dir path; directory containing only the FASTQ sRNAseq 
#' samples for analysis. Note that all samples in this directory will be used by 
#' this function. 
#' @param condaenv character; name or directory of the Conda environment to use
#' where `ShortStack` (>4.0) is installed
#' 
#' @param output_dir path; directory to store output. 
#' @param genomefile path; path to a FASTA file. 
#' @param threads numeric; set the number of threads to use where more threads 
#' means a faster completion time. Default is 6. 
#' 
#' @param mmap character; define how to handle multi-mapped reads. Choose from 
#' "u", "f" or "r". Where "u" means only uniquely-aligned reads are used as 
#' weights for placement of multi-mapped reads. Where "f" means fractional 
#' weighting scheme for placement of multi-mapped reads and "r" mean 
#' multi-mapped read placement is random. Default is "u". 
#' 
#' @param dicermin integer; the minimum size in nucleotides of a valid small 
#' RNA. This option sets the bounds to discriminate dicer-derived small RNA loci 
#' from other loci. Default is 20.
#' @param dicermax integer; the maximum size in nucleotides of a valid small 
#' RNA. This option sets the bounds to discriminate dicer-derived small RNA loci 
#' from other loci.  Default is 24.
#' 
#' @param mincov numeric; minimum alignment depth, in units of reads per 
#' million, required to nucleate a small RNA cluster during de novo cluster 
#' search. Must be an floating point number > 0. Default is 2.
#' 
#' @param pad integer; initial peaks are merged if they are this distance or 
#' less from each other. Must >= 1, default is 75. 
#' @param tidy logical; removes unnecessary extra output files when set to TRUE. 
#' 
#' 
#' @return
#' The OS commands generate output into the users desired location, generating
#' two folders:
#' 
#' * "1_de_novo_detection" : Stores output from the detection of de novo sRNA-producing loci
#' * "2_alignment_results" : Stores output from alignment 
#' 
#' 
#' The first folder stores the novo detection of sRNA-producing loci and
#' alignment, where the output of each are stored in their respective folders in 
#' the users desired location. The de novo detection of sRNA-producing loci
#' analyses each sample to identify de novo sRNA-producing loci (ie. sRNA
#' clusters), and joins these results into a single file called "locifile.txt". 
#' The second folder stores the alignment and clustering results for each sample 
#' to the genome reference along with the file containing the de novo sRNA 
#' clusters. The output of each sample is stored within its own sample folder. 
#' 
#' The alignment results generate one folder per sample, where the results 
#' are stored. As default this includes the alignment files (BAM) and 
#' sRNA results (.txt). The sRNA results (.txt) are imported into R for 
#' downstream analysis by utilizing the [mobileRNA::RNAimport()] function. 
#' 
#' 
#' The function generates a number of extra files for each sample and are not 
#' required for the downstream analysis. See `ShortStack` documentation for more
#' information (\link{https://github.com/MikeAxtell/ShortStack}). As default
#' these files are deleted. This is determined by the `tidy` argument. 
#' 
#' @examples \dontrun{
#' 
#' samples <- file.path(system.file("extdata",package="mobileRNA"))
#' 
#' output_location <- tempdir()
#' 
#' sRNAmapper(input_files_dir = samples, 
#' output_dir = output_location, 
#' genomefile = output_assembly_file,
#' condaenv = "ShortStack4")
#' }
#' 
#' @importFrom reticulate "use_condaenv" 
#' @importFrom GenomicRanges GRangesList
#' @importFrom Repitools annoGR2DF
#' @importFrom rtracklayer import.gff
#' @importFrom GenomicRanges reduce
#' @importFrom utils write.table
#' 
#' @export
sRNAmapper <- function(input_files_dir, output_dir, 
                       genomefile, 
                       condaenv,
                       threads = 6, 
                       mmap = "u", 
                       dicermin = 20,
                       dicermax = 24, 
                       mincov = 2, 
                       pad = 75, 
                       tidy = TRUE){
  
  # 1 - Verify that the OS is either Linux or macOS
  os <- tolower(Sys.info()[['sysname']])
  if (os != "linux" && os != "darwin") {
    stop("ShortStack can only run on Linux or macOS.")
  }
  # set conda envrioment of shortstack
  reticulate::use_condaenv(condaenv, required = TRUE)
  exists_res <- shortstack_exists()
  if(exists_res < 4 | is.null(exists_res)){
    stop("ShortStack application :
    --- ShortStack is either not installed or installed incorrectly
    --- Or the version is too old, please ensure most updated version is installed.")
  }
  # 3 - generate output folders (check if they exist )
  path_1 <- file.path(output_dir, "1_de_novo_detection")
  path_2 <- file.path(output_dir, "2_alignment_results")
  if (!dir.exists(path_1)) {
    dir.create(path_1)
  }
  if (!dir.exists(path_2)) {
    dir.create(path_2)
  }
 
  # run command
  names_input_files_dir <- list.files(input_files_dir)
  message("mapping with ShortStack ...")
  stats <- file.path(path_1, "stats_alignment.txt")
  system(paste0(">> ", stats))
 
  for (i in seq_along(names_input_files_dir)){
    # file name:
    readfile_name <- sub("\\.\\w+$", "", names_input_files_dir[i])
    # step 1 - map as per 
    # file out put: 
    file_outdir <- file.path(path_1, readfile_name)
    
    shortstack_cmd_1 <- c(
      shQuote(Sys.which("shortstack")),
      "--readfile", shQuote(file.path(input_files_dir,names_input_files_dir[i])),
      "--genomefile", shQuote(genomefile), 
      "--threads", shQuote(threads),
      "--mmap", shQuote(mmap),
      "--dicermin", shQuote(dicermin),
      "--dicermax", shQuote(dicermax),
      "--nohp", 
      "--mincov",shQuote(mincov),
      "--pad", shQuote(pad), 
      "--outdir", shQuote(file_outdir),">>", shQuote(stats), "2>&1"
    ) 
    
    shortstack_cmd_1 <- paste(shortstack_cmd_1,collapse = " ")
    shortstack_cmd_1 <- gsub("^ *| *$", "", shortstack_cmd_1)
    
    # Run ShortStack using system command
    system(shortstack_cmd_1, intern=FALSE)
  }
 
  # 5 - merge loci into file
  map_1_files_loci <- list.dirs(path_1, full.names = TRUE, recursive = TRUE)
  map_1_files_loci <- map_1_files_loci[!map_1_files_loci == path_1]
  samples <- basename(map_1_files_loci)
  
  gff_alignment <- GenomicRanges::GRangesList()
  for (i in samples) {
    file_path <- file.path(path_1, i, "Results.gff3")
    if (file.exists(file_path)) {
    gff_alignment[[i]] <- rtracklayer::import.gff(file_path)
    } else{
      Stop("File does not exist:", file_path, "\n")
    }
  }
  
  gff_merged <- GenomicRanges::reduce(unlist(gff_alignment), 
                                      ignore.strand = TRUE)
  gff_merged <- Repitools::annoGR2DF(gff_merged)
  locifile_txt <- data.frame(Locus = paste0(gff_merged$chr, ":", 
                                            gff_merged$start,"-", 
                                            gff_merged$end), 
                             Cluster = paste0("cluster_", 
                                              seq_len(nrow(gff_merged))))
  
  loci_out <- file.path(path_1,"locifile.txt")
  utils::write.table(locifile_txt, file = loci_out, quote = FALSE, 
                     sep = "\t", row.names = FALSE, col.names = TRUE)
   
  # 6 - mapping 2
  stats_2 <- file.path(path_2, "stats_alignment.txt")
  system(paste0(">> ", stats_2))
  
  for (i in seq_along(names_input_files_dir)){
    # file name:
    readfile_name <- sub("\\.\\w+$", "", names_input_files_dir[i])
    # step 1 - map as per 
    # file out put: 
    file_outdir <- file.path(path_2, readfile_name)
    
    shortstack_cmd_2 <- c(
      shQuote(Sys.which("shortstack")),
      "--readfile", shQuote(file.path(input_files_dir,names_input_files_dir[i])),
      "--locifile", shQuote(loci_out), 
      "--genomefile", shQuote(genomefile), 
      "--threads", shQuote(threads),
      "--mmap", shQuote(mmap),
      "--dicermin", shQuote(dicermin),
      "--dicermax", shQuote(dicermax),
      "--nohp", 
      "--mincov",shQuote(mincov),
      "--pad", shQuote(pad), 
      "--outdir", shQuote(file_outdir),">>", shQuote(stats_2), "2>&1"
    ) 
    
    shortstack_cmd_2 <- paste(shortstack_cmd_2,collapse = " ")
    shortstack_cmd_2 <- gsub("^ *| *$", "", shortstack_cmd_2)
    
    # Run ShortStack using system command
    system(shortstack_cmd_2, intern=FALSE)
  }
 
    # remove excess files 
  if (tidy) {
    # List all directories in path_1 and exclude path_1 itself
    map_1_files <- list.dirs(path_1, full.names = TRUE, recursive = TRUE)
    map_1_files <- map_1_files[!map_1_files == path_1]
    
    # Iterate over each directory and remove files that don't match the conditions
    for (j in seq_along(map_1_files)) {
      rm_cmd_1 <- paste0("find '", map_1_files[j], "' -type f ! -name '*.bam' -exec rm {} \\;")
      system(rm_cmd_1, intern=FALSE)
    }
    
    # List all directories in path_2 and exclude path_2 itself
    map_2_files <- list.dirs(path_2, full.names = TRUE, recursive = TRUE)
    map_2_files <- map_2_files[!map_2_files == path_2]
    
    # Iterate over each directory and remove files that don't match the conditions
    for (k in seq_along(map_2_files)) {
      rm_cmd_2 <- paste0("find '", map_2_files[k], "' -type f ! \\( -name '*.bam' -o -name 'Results.txt' \\) -exec rm {} \\;")
      system(rm_cmd_2, intern=FALSE)
    }
  }
 
  cat("\n")
  cat("\n")
  message(" --- Mapping of sRNAseq samples is complete --- ")
  message("Results saved to: ", path_1, " & ", path_2)
  message("Loci file saved to: ", loci_out)
  
}
