#' Map & cluster sRNAseq reads using `ShortStack`
#'
#' @description Map sRNAseq FASTQ samples to a given genome assembly using 
#' ShortStack. 
#' 
#' @details
#' This function undertakes our recommended workflow for mapping sRNAseq reads
#' to a given genome for the analysis of small RNAs in a biological system. 
#' 
#' The function invokes a number of OS commands, and is dependent on the 
#' installation of `ShortStack` and it's dependencies within the same 
#' environment. Please refer to \link{} for more information. Please not that 
#' `ShortStack` is only compartible with Linux and Mac operating systems. 
#' 
#' It uses `ShortStack` to map and cluster sRNAs for each sample 
#' within a directory. Within the output directory (`output_dir`), the function
#' generates two new directories called  "1_map" and "2_map", plus, a plain
#' text file called "locifile.txt". The function undertakes an initial mapping 
#' step in order to generate a list of de novo sRNA-producing loci. The results 
#' for each sample are saved to "1_map", producing one folder for each sample. 
#' The de novo sRNA-producing loci identified in each sample are merged; 
#' generating the "locifile.txt" file. The second mapping step utilises the 
#' "locifile.txt" file to optimise the analysis. The results of mapping step two
#' are saved to the "2_map" directory, where one folder for each sample is 
#' generated and stores the results. 
#' 
#' The function generates a number of extra files for each sample and are not 
#' required for the downstream analysis. Hence, as default these files are 
#' deleted. This is determined by the `tidy` argument. 
#' 
#'
#' @param input_files_dir path; directory containing only the FASTQ sRNAseq 
#' samples for analysis. Note that all samples in this directory will be used by 
#' this function. 
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
#' from other loci.  Default is 25.
#' 
#' @param mincov numeric; minimum alignment depth, in units of reads per 
#' million, required to nucleate a small RNA cluster during de novo cluster 
#' search. Must be an floating point number > 0. Default is 2.
#' 
#' @param pad interger; initial peaks are merged if they are this distance or 
#' less from each other. Must >= 1, default is 75. 
#' @param tidy logical; removes unnecessary extra output files when set to TRUE. 
#' 
#' @examples \dontrun{
#' 
#' options(shortstack_path = "[path to shortstack]")
#' 
#' samples <- file.path(system.file("extdata",package="mobileRNA"))
#' 
#' output_location <- tempdir()
#' 
#' sRNAmapper(input_files_dir = samples, 
#' output_dir = output_location, 
#' genomefile = output_assembly_file)
#' }
#' 
#' 
#' 
#'
sRNAmapper <- function(input_files_dir, output_dir, 
                       genomefile,
                       threads = 6, 
                       mmap = "u", 
                       dicermin = 20,
                       dicermax = 25, 
                       mincov = 2, 
                       pad = 75, 
                       tidy = TRUE){
  
  # 1 - Verify that the OS is either Linux or macOS
  os <- tolower(Sys.info()[['sysname']])
  if (os != "linux" && os != "darwin") {
    stop("ShortStack can only run on Linux or macOS.")
  }
  
  system_dependencies <- c("samtools", "ShortStack", "bowtie", "gzip")
  results <- sapply(system_dependencies, function(dep) {
    cmd <- paste("conda list | grep ", dep)
    result <- system(cmd, intern =TRUE)
    # Check if the result contains the dependency name
    return(grepl(dep, result))
  })
  names(results) <- system_dependencies
  if(length(unlist(results)) != 4){
    stop("Not all system dependancies found. Please create a environment holding 
         ShortStack and it's dependancies; samtools, bowtie and gzip.")
  }
  
  # 3 - generate output folders (check if they exist )
  path_1 <- file.path(output_dir, "1_map")
  path_2 <- file.path(output_dir, "2_map")
  if (!dir.exists(path_1)) {
    dir.create(path_1)
  }
  if (!dir.exists(path_2)) {
    dir.create(path_2)
  }
  
  # 4 - mapping 1 
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
      shQuote(getOption("shortstack_path")),
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
    system(shortstack_cmd_1, intern=FALSE, show.output.on.console=FALSE)
  }
  message("Completed mapping step 1. Results saved to: ", path_1)
  
  # 5 - merge loci into file
  map_1_files_loci <- list.dirs(path_1, full.names = TRUE, recursive = TRUE)
  map_1_files_loci <- map_1_files_loci[!map_1_files_loci == path_1]
  samples <- basename(map_1_files_loci)
  
  gff_alignment <- GenomicRanges::GRangesList()
  for (i in samples) {
    gff_alignment[[i]] <- rtracklayer::import.gff(file.path(path_1, samples[i], 
                                                         "ShortStack_All.gff3"))
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
  
  message("Saved sRNA loci to: ", loci_out)
  
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
      shQuote(getOption("shortstack_path")),
      "--readfile", shQuote(file.path(input_files_dir,names_input_files_dir[i])),
      "--locifile", shQuote(locifile_txt), 
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
    system(shortstack_cmd_2, intern=FALSE, show.output.on.console=FALSE)
  }
  message("Completed mapping step 2. Results saved to: ", path_2)
  
  # remove excess files 
  if (tidy) {
    # List all directories in path_1 and exclude path_1 itself
    map_1_files <- list.dirs(path_1, full.names = TRUE, recursive = TRUE)
    map_1_files <- map_1_files[!map_1_files == path_1]
    
    # Iterate over each directory and remove files that don't match the conditions
    for (j in seq_along(map_1_files)) {
      rm_cmd_1 <- paste0("find '", map_1_files[j], "' -type f ! -name '*.bam' -exec rm {} \\;")
      system(rm_cmd_1, intern=FALSE, show.output.on.console=FALSE)
    }
    
    # List all directories in path_2 and exclude path_2 itself
    map_2_files <- list.dirs(path_2, full.names = TRUE, recursive = TRUE)
    map_2_files <- map_2_files[!map_2_files == path_2]
    
    # Iterate over each directory and remove files that don't match the conditions
    for (k in seq_along(map_2_files)) {
      rm_cmd_2 <- paste0("find '", map_2_files[k], "' -type f ! \\( -name '*.bam' -o -name 'Results.txt' \\) -exec rm {} \\;")
      system(rm_cmd_2, intern=FALSE, show.output.on.console=FALSE)
    }
  }
  return("Mapping of sRNAseq samples is complete.")
}
