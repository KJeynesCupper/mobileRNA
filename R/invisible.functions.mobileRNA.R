#------------------------------------------------------------#
# Title:  Invisible functions                                #
# Author: Katie Jeynes-Cupper (kej031@student.bham.ac.uk) #
# Date:   01.02.23                                           #
#------------------------------------------------------------#
################ remove mapping errors (RNAmobile function) ####################
.remove_mapping_errors <- function(data, controls) {
  class_colnames  <- data %>% dplyr::select(paste0("Count_", controls))
  x <- c()
  if (length(colnames(class_colnames)) > 1){
    for (j in seq_len(nrow(data))){
      if(rowSums(class_colnames[j,])>0){
        x <- c(x,j)
      }else 
        x <- c(x)
    }
  } else
    if(length(colnames(class_colnames)) == 1){
      x <- c()
      for (k in seq_len(nrow(data)) ){
        if(stats::na.omit(as.numeric(data[k,colnames(class_colnames)],
                                     na.rm=TRUE))!= 0){
          x <- c(x,k)
        }
      }
    }
  if(is.null(x)){
    data <- data
  } else
    data <- data[-x,]
  return(data)
}
################ Remove mapping errors  #########################
.remove_mapping_errors_V2 <- function(data,  controls, genome.ID) {
    if (base::missing(controls) || !base::inherits(controls, "character")) {
    stop("Please specify a character vector storing names of control 
          replicates")
  }
    if (base::missing(genome.ID) || genome.ID %in% "") {
    stop("Please specify a single character string which is present in the all 
           the chromosomes within the foriegn genome")
  }
    data_native <- data %>% dplyr::filter(!grepl(genome.ID,chr))
  # subset data to find all rows of forign genome
    data_select <- data %>% dplyr::filter(grepl(genome.ID,chr))
    class_colnames  <- data_select %>% dplyr::select(paste0("Count_", controls))
    if (length(colnames(class_colnames)) > 1){
    x <- c()
    for (j in seq_len(nrow(data_select)) ){
    if(sum(stats::na.omit(as.numeric(data_select[j,colnames(class_colnames)],
                                        na.rm=TRUE)))>0){
    x <- c(x,j)
      }
    }
  } else
    if (length(colnames(class_colnames)) == 1){
    x <- c()
    for (k in seq_len(nrow(data_select)) ){
    if(stats::na.omit(as.numeric(data_select[k,colnames(class_colnames)],
                                     na.rm=TRUE))!= 0){
          x <- c(x,k)
        }
      }
    }
    if(is.null(x)){
    data_id <- data_select
  } else
    if(!is.null(x)){ 
      data_id <- data_select[-x,]
    }
    data <- rbind(data_native,data_id)
}

################ DESE2 function (RNAdifferentialAnalysis function) #############
.DESeq_normalise <- function(data, conditions){
  column.data <- data.frame(conditions=as.factor(conditions))
  base::rownames(column.data) <- base::colnames(data)
  count.data.set <- DESeq2::DESeqDataSetFromMatrix(countData=data,
                                                   colData=column.data,
                                                   design= ~ conditions)
  count.data.set$conditions <- stats::relevel(count.data.set$conditions,
                                              conditions[1])
  out <- DESeq2::estimateSizeFactors(count.data.set)
  return(out)
}
################ EDGER function (RNAdifferentialAnalysis function) #############
.edgeR_normalise <- function(data, conditions){
  d <- edgeR::DGEList(counts = data, group = factor(conditions))
  result <- edgeR::calcNormFactors(d)
  result$samples
  result<- edgeR::estimateDisp(result)
  result$common.dispersion
  return(result)
}

################ Find RNA complementary sequence (RNAsequence function) ########
find_complementary_sequenceRNA <- function(seq) {
  # conversions
  conversion_nucleotides <- c(A = "U", U = "A", C = "G", G = "C")

  # calculate complementary nt for each
  complementary_calc <- sapply(strsplit(seq, ""), function(nucleotide) {
    conversion_nucleotides[nucleotide]
  })

  # Combine into string
  output <- paste0(complementary_calc, collapse = "")
  return(output)
}

################ Find DNA complementary sequence (RNAsequence function) #######
find_complementary_sequenceDNA <- function(seq) {
  # conversions
  conversion_nucleotides <- c(A = "T", U = "A", C = "G", G = "C")

  # calc complementary nt for each in string
  complementary_calc <- sapply(strsplit(seq, ""), function(nucleotide) {
    conversion_nucleotides[nucleotide]
  })

  # Combine into string
  output <- paste0(complementary_calc, collapse = "")
  return(output)
}

##################shortstack exists #####################################

shortstack_exists <- function(){
  # Replace "your_package_name" with the name of the package you want to check
  package_name <- "shortstack"
  cmd <- c(shQuote(package_name),  "--version")
  cmd_a <- paste(cmd,collapse = " ")
  cmd_b <- gsub("^ *| *$", "", cmd_a)
  out <- system(cmd_b, intern=TRUE )
  version <- gsub(".*?([0-9.]+).*", "\\1", out)
  return(version)
}

################## import gff #####################################

gff_import <- function(gff_file, nrows = -1) {
  # Read the GFF file into a dataframe
  gff <- utils::read.table(
    gff_file,
    sep = "\t",
    as.is = TRUE,
    quote = "",
    header = FALSE,
    comment.char = "#",
    nrows = nrows,
    colClasses = c("character", "character", "factor", "integer",
                   "integer", "character", "character", "character","character")
  )
  # convert problematic symbol error
  gff$V9 <- gsub('%', '=',gff$V9)
  colnames(gff) <- c("seqname", "source", "type", "start", "end",
                     "score", "strand", "phase", "attributes")
  
  if (any(is.na(gff$start)) || any(is.na(gff$end))) {
    stop("The 'start' or 'end' columns contain missing values.")
  }
  return(gff)
}



################### convert character to factor in gramges #####

convertChar2Factor <- function(gr) {
  charCols <- sapply(S4Vectors::elementMetadata(gr), is.character)
  
  if (any(charCols)) {
    gr_metadata <- S4Vectors::elementMetadata(gr)
    gr_metadata[charCols] <- lapply(gr_metadata[charCols], as.factor)
    S4Vectors::elementMetadata(gr) <- gr_metadata
  }
  metadata_cols <- S4Vectors::elementMetadata(gr)
  charlist_cols <- sapply(metadata_cols@listData, function(col) 
    inherits(col, "CompressedCharacterList"))
  
 if (any(charlist_cols)) {
   metadata_cols <- metadata_cols[!charlist_cols]
   S4Vectors::elementMetadata(gr) <- metadata_cols
  }
  return(gr)
}






################## core map #####################################
core_map <- function(input_files_dir, output_dir, genomefile, condaenv, 
             threads,mmap, dicermin,dicermax, mincov, pad, tidy){
  
  
  # 3 - generate output folders (check if they exist )
  path_1 <- file.path(output_dir, "1_de_novo_detection")
  path_2 <- file.path(output_dir, "2_sRNA_results")
  if (!dir.exists(path_1)) {
    dir.create(path_1)
  }
  if (!dir.exists(path_2)) {
    dir.create(path_2)
  }
  
  # run command
  names_input_files_dir <- list.files(input_files_dir)
  message("mapping with ShortStack ...")
  stats <- file.path(path_1, "log.txt")
  system(paste0(">> ", stats))
  
  # 1 - bowtie - build 
  # check is alread y done
  index_extension <- "ebwt"
  base_genomefile <- tools::file_path_sans_ext(basename(genomefile))
  dir_genomfile <- dirname(genomefile)
  files_dir <- c("ls", dir_genomfile)
  files_dir <- paste(files_dir,collapse = " ")
  files_dir <- gsub("^ *| *$", "", files_dir)
  files_dir_res <- system(files_dir, intern=TRUE)
  
  index_found <- grep(paste("^", base_genomefile, ".*\\.", index_extension, 
                            "$", sep = ""), files_dir_res)
  base_genomefile_path <- tools::file_path_sans_ext(genomefile)
  
  if (length(index_found) > 0) {
    message("Genome index already built.")
    nline <- paste("echo Genome index already built. >> ",  shQuote(stats))
    system(nline, intern=FALSE)
    nline <- paste("echo '' >> ",  shQuote(stats))
    system(nline, intern=FALSE)
  } else {
    message("Genome index not already built. Building with Bowtie.")
    nline <- paste("echo Genome index not already built. Building with Bowtie: >> ",  shQuote(stats))
    system(nline, intern=FALSE)
    nline <- paste("echo '' >> ",  shQuote(stats))
    system(nline, intern=FALSE)
    
    bowtie_build_cmd <- c("bowtie-build --threads", threads, 
                          genomefile, base_genomefile_path, ">>", 
                          shQuote(stats), "2>&1") 
    bowtie_build_cmd <- paste(bowtie_build_cmd,collapse = " ")
    bowtie_build_cmd <- gsub("^ *| *$", "", bowtie_build_cmd)
    system(bowtie_build_cmd, intern=FALSE) ## -- get it to print to summarydoc. 
    cat("\n")
    message("Genome index build. ")
    nline <- paste("echo Genome index complete! >> ",  shQuote(stats))
    system(nline, intern=FALSE)
    nline <- paste("echo '' >> ",  shQuote(stats))
    system(nline, intern=FALSE)
  }
  
  for (i in seq_along(names_input_files_dir)){
    # file name:
    readfile_name <- sub("\\.\\w+$", "", names_input_files_dir[i])
    # step 1 - map as per 
    # file out put: 
    file_outdir <- file.path(path_1, readfile_name)
    cat("\n", file = stats, append = TRUE)
    time_cmd <- paste(c("echo", shQuote(as.character(Sys.time())), 
                        "Working with sample:", 
                        shQuote(readfile_name),
                        ">>", shQuote(stats)), collapse = " ")
    time_cmd <- gsub("^ *| *$", "", time_cmd)
    system(time_cmd, intern=FALSE)
    
    shortstack_cmd_1 <- c(
      shQuote(Sys.which("shortstack")),
      "--readfile", shQuote(file.path(input_files_dir,
                                      names_input_files_dir[i])),
      "--genomefile", shQuote(genomefile), 
      "--threads", shQuote(threads),
      "--mmap", shQuote(mmap),
      "--dicermin", shQuote(dicermin),
      "--dicermax", shQuote(dicermax),
      "--nohp", 
      "--mincov 0.5",
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
      stop("File does not exist:", file_path, "\n")
    }
  }
  gff_merged <- GenomicRanges::reduce(unlist(gff_alignment), 
                                      ignore.strand = TRUE)
  gff_merged <- as.data.frame(gff_merged)
  colnames(gff_merged)[1] <- "chr"
  if('*' %in% gff_merged$strand){
    gff_merged <- gff_merged[, -match("strand", colnames(gff_merged))]
  }
  
  locifile_txt <- data.frame(Locus = paste0(gff_merged$chr, ":", 
                                            gff_merged$start,"-", 
                                            gff_merged$end), 
                             Cluster = paste0("cluster_", 
                                              seq_len(nrow(gff_merged))))
  
  loci_out <- file.path(path_1,"locifile.txt")
  utils::write.table(locifile_txt, file = loci_out, quote = FALSE, 
                     sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # 6 - mapping 2
  stats_2 <- file.path(path_2, "log.txt")
  system(paste0(">> ", stats_2))
  ###################.  check bam location ull is correct 
  for (i in seq_along(names_input_files_dir)){
    # file name:
    readfile_name <- sub("\\.\\w+$", "", names_input_files_dir[i])
    # step 1 - map as per 
    # file out put: 
    file_outdir <- file.path(path_2, readfile_name)
    cat("\n", file = stats_2, append = TRUE)
    time_cmd <- paste(c("echo", shQuote(as.character(Sys.time())), 
                        "Working with sample:", 
                        shQuote(readfile_name),
                        ">>", shQuote(stats_2)), collapse = " ")
    time_cmd <- gsub("^ *| *$", "", time_cmd)
    system(time_cmd, intern=FALSE)
    shortstack_cmd_2 <- c(
      shQuote(Sys.which("shortstack")),
      "--bamfile", shQuote(file.path(path_1,readfile_name,
                                     paste0(readfile_name,".bam"))),
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
    
    # Iterate over each directory and remove files that don't match the 
    # conditions
    for (j in seq_along(map_1_files)) {
      rm_cmd_1 <- paste0("find '", map_1_files[j],
                         "' -type f ! -name '*.bam' -exec rm {} \\;")
      system(rm_cmd_1, intern=FALSE)
    }
    
    # List all directories in path_2 and exclude path_2 itself
    map_2_files <- list.dirs(path_2, full.names = TRUE, recursive = TRUE)
    map_2_files <- map_2_files[!map_2_files == path_2]
    
    # Iterate over each directory and remove files that don't match the 
    # conditions
    for (k in seq_along(map_2_files)) {
      rm_cmd_2 <- paste0("find '", map_2_files[k], 
    "' -type f ! \\( -name '*.bam' -o -name 'Results.txt' \\) -exec rm {} \\;")
      system(rm_cmd_2, intern=FALSE)
    }
  }
  
  cat("\n")
  cat("\n")
  message(" --- Mapping of sRNAseq samples is complete --- ")
  message("Results saved to: ", path_1, " & ", path_2)
  message("Loci file saved to: ", loci_out)
  
}



################## mobile map  #####################################
mobile_map <- function(input_files_dir, output_dir, genomefile, condaenv, 
                     threads,mmap, dicermin,dicermax, mincov, pad, tidy){
  
  # 3 - generate output folders (check if they exist )
  path_1 <- file.path(output_dir, "1_de_novo_detection")
  path_2 <- file.path(output_dir, "2_sRNA_results")
  if (!dir.exists(path_1)) {
    dir.create(path_1)
  }
  if (!dir.exists(path_2)) {
    dir.create(path_2)
  }
  message("Begining pre-processing for mobile sRNAseq ...")
  
  stats <- file.path(path_1, "log.txt")
  system(paste0(">> ", stats))
  
  # 1 - bowtie - build 
  # check is alread y done
  index_extension <- "ebwt"
  base_genomefile <- tools::file_path_sans_ext(basename(genomefile))
  dir_genomfile <- dirname(genomefile)
  files_dir <- c("ls", dir_genomfile)
  files_dir <- paste(files_dir,collapse = " ")
  files_dir <- gsub("^ *| *$", "", files_dir)
  files_dir_res <- system(files_dir, intern=TRUE)
  
   index_found <- grep(paste("^", base_genomefile, ".*\\.", index_extension, 
                             "$", sep = ""), files_dir_res)
   base_genomefile_path <- tools::file_path_sans_ext(genomefile)

  if (length(index_found) > 0) {
    message("Genome index already built.")
    nline <- paste("echo Genome index already built. >> ",  shQuote(stats))
    system(nline, intern=FALSE)
    nline <- paste("echo '' >> ",  shQuote(stats))
    system(nline, intern=FALSE)
  } else {
    message("Genome index not already built. Building with Bowtie.")
    nline <- paste("echo Genome index not already built. Building with Bowtie: >> ",  shQuote(stats))
    system(nline, intern=FALSE)
    nline <- paste("echo '' >> ",  shQuote(stats))
    system(nline, intern=FALSE)
    
    bowtie_build_cmd <- c("bowtie-build --threads", threads, 
                          genomefile, base_genomefile_path, ">>", 
                          shQuote(stats), "2>&1") 
    bowtie_build_cmd <- paste(bowtie_build_cmd,collapse = " ")
    bowtie_build_cmd <- gsub("^ *| *$", "", bowtie_build_cmd)
    system(bowtie_build_cmd, intern=FALSE) ## -- get it to print to summarydoc. 
    cat("\n")
    message("Genome index build. ")
    nline <- paste("echo Genome index complete! >> ",  shQuote(stats))
    system(nline, intern=FALSE)
    nline <- paste("echo '' >> ",  shQuote(stats))
    system(nline, intern=FALSE)
  }
  
  # 2-  run bowtie mappling ## -- get it to print to a summary doc. 
  names_input_files_dir <- list.files(input_files_dir)
  message("mapping with Bowtie ...")
  
  p1_path <- file.path(path_1,"Alignment")
  p1_dir <- dir.create(p1_path)
  
  #bw message 
  bowtie_cmd_message <-c("Mapping with Bowtie, Running Command: 
  bowtie -p", threads, 
  "-v 0 -a -m 1 --sam -x [genomefile] [fastqfile] | samtools view -bS",
  "| samtools sort -o [outputfile]")
  bowtie_cmd_message <- paste(bowtie_cmd_message,collapse = " ")
  bowtie_cmd_message <- gsub("^ *| *$", "", bowtie_cmd_message)
  cat(bowtie_cmd_message, file = stats, append = TRUE)
  cat("\n", file = stats, append = TRUE)
  
  for (i in seq_along(names_input_files_dir)){
    # file name:
    cat("\n", file = stats, append = TRUE)
    readfile_name <- sub("\\.\\w+$", "", names_input_files_dir[i])
    # step 1 - map as per 
    # file output folder: 
    op <- paste0(readfile_name,".bam")
    output_file <- file.path(path_1,"Alignment",op)
    time_cmd <- paste(c("echo", shQuote(as.character(Sys.time())), 
                        "Alignment with Bowtie for sample:", 
                          shQuote(readfile_name),  ">>", 
                        shQuote(stats)), collapse = " ")
    time_cmd <- gsub("^ *| *$", "", time_cmd)
    system(time_cmd, intern=FALSE)
    
    bowtie_cmd_1 <- c("( bowtie",
                      "-p", shQuote(threads), 
                      "-v 0 -a -m 1 --sam -x", 
                      shQuote(base_genomefile_path),
                      shQuote(file.path(input_files_dir,
                                        names_input_files_dir[i])), 
                      "| samtools view -bS", 
                      "| samtools sort -o",shQuote(output_file), ")",
                      "2>>", shQuote(stats))
    
    bowtie_cmd_1 <- paste(bowtie_cmd_1,collapse = " ")
    bowtie_cmd_1 <- gsub("^ *| *$", "", bowtie_cmd_1)
    system(bowtie_cmd_1, intern=FALSE)
    
  }
  
  #bw message 
  cat("\n", file = stats, append = TRUE)
  clustering_cmd_message <-c("echo Running sRNA De Novo Detection: >>", 
                             shQuote(stats))
  
  clustering_cmd_message <- paste(clustering_cmd_message,collapse = " ")
  clustering_cmd_message <- gsub("^ *| *$", "", clustering_cmd_message)
  system(clustering_cmd_message, intern=FALSE)
  
  
  mapped_file_names <- list.files(p1_path)
  p2_path <- file.path(path_1,"Cluster")
  p2_dir <- dir.create(p2_path)

  # 3 - shorstack clustering 
  for (i in seq_along(mapped_file_names)){
    # file name:
    readfile_name <- sub("\\.\\w+$", "", mapped_file_names[i])
    # file out put: 
    file_outdir <- file.path(p2_path, readfile_name)
    cat("\n", file = stats, append = TRUE)
    time_cmd <- paste(c("echo", shQuote(as.character(Sys.time())), 
                        "De novo sRNA detection for sample:", 
                        shQuote(readfile_name),
                        ">>", shQuote(stats)), collapse = " ")
    time_cmd <- gsub("^ *| *$", "", time_cmd)
    system(time_cmd, intern=FALSE)
  
    shortstack_cmd_1 <- c(
      shQuote(Sys.which("shortstack")),
      "--bamfile", shQuote(file.path(p1_path,mapped_file_names[i])),
      "--genomefile", shQuote(genomefile), 
      "--threads", shQuote(threads),
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
  
  # 5 - merge loci into file. clustered_files
  clustered_files <- list.dirs(p2_path, full.names = TRUE, recursive = TRUE)
  clustered_files <- clustered_files[!clustered_files == p2_path]
  samples <- basename(clustered_files)
  
  gff_alignment <- GenomicRanges::GRangesList()
  for (i in samples) {
    file_path <- file.path(p2_path, i, "Results.gff3")
    if (file.exists(file_path)) {
      gff_alignment[[i]] <- rtracklayer::import.gff(file_path)
    } else{
      stop("File does not exist:", file_path, "\n")
    }
  }
  gff_merged <- GenomicRanges::reduce(unlist(gff_alignment), 
                                      ignore.strand = TRUE)
  gff_merged <- as.data.frame(gff_merged)
  colnames(gff_merged)[1] <- "chr"
  if('*' %in% gff_merged$strand){
    gff_merged <- gff_merged[, -match("strand", colnames(gff_merged))]
  }
  locifile_txt <- data.frame(Locus = paste0(gff_merged$chr, ":", 
                                            gff_merged$start,"-", 
                                            gff_merged$end), 
                             Cluster = paste0("cluster_", 
                                              seq_len(nrow(gff_merged))))
  
  loci_out <- file.path(path_1,"locifile.txt")
  utils::write.table(locifile_txt, file = loci_out, quote = FALSE, 
                     sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # move bam files into respecitve dir. 
  for(i in seq_along(mapped_file_names)){
    bam_path <- file.path(p1_path,mapped_file_names[i])
    file_name_without_ext <- tools::file_path_sans_ext(mapped_file_names[i])
    saveloc <- file.path(p2_path, file_name_without_ext)
    
    organise_cmd <- c("mv", bam_path, paste0(saveloc, "/"))
    organise_cmd <- paste(organise_cmd,collapse = " ")
    organise_cmd <- gsub("^ *| *$", "", organise_cmd)
    system(organise_cmd, intern=FALSE)
  }
  
  ###### move and remove cluster dir. 
  move_clusters_cmd <- c("mv", paste0(p2_path, "/*"), paste0(path_1, "/"))
  move_clusters_cmd <- paste(move_clusters_cmd,collapse = " ")
  move_clusters_cmd <- gsub("^ *| *$", "", move_clusters_cmd)
  system(move_clusters_cmd, intern=FALSE)
  
  # delete the alignment & cluster folder folder:
  delete_align <- c("rm -r", p1_path, p2_path)
  delete_align <- paste(delete_align,collapse = " ")
  delete_align <- gsub("^ *| *$", "", delete_align)
  system(delete_align, intern=FALSE)
  
  ####### final clustering 
  stats_2 <- file.path(path_2, "log.txt")
  system(paste0(">> ", stats_2))
  
  clustering_cmd_message<-c("echo Running sRNA Clustering with loci information >>", 
                             shQuote(stats_2))
  
  clustering_cmd_message <- paste(clustering_cmd_message,collapse = " ")
  clustering_cmd_message <- gsub("^ *| *$", "", clustering_cmd_message)
  system(clustering_cmd_message, intern=FALSE)
  
  for (i in seq_along(mapped_file_names)){
    # file name:
    readfile_name <- sub("\\.\\w+$", "", mapped_file_names[i])
    # step 1 - map as per 
    # file out put: 
    file_outdir <- file.path(path_2, readfile_name)
    cat("\n", file = stats_2, append = TRUE)
    time_cmd <- paste(c("echo", shQuote(as.character(Sys.time())), 
                        "Working on sample:", 
                        shQuote(readfile_name),
                        ">>", shQuote(stats_2)), collapse = " ")
    time_cmd <- gsub("^ *| *$", "", time_cmd)
    system(time_cmd, intern=FALSE)
    
    shortstack_cmd_2 <- c(
      shQuote(Sys.which("shortstack")),
      "--bamfile", shQuote(file.path(path_1,readfile_name, 
                                     mapped_file_names[i])),
      "--locifile", shQuote(loci_out), 
      "--genomefile", shQuote(genomefile), 
      "--threads", shQuote(threads),
      "--dicermin", shQuote(dicermin),
      "--dicermax", shQuote(dicermax),
      "--mincov",shQuote(mincov),
      "--pad", shQuote(pad), 
      "--outdir", shQuote(file_outdir),">>", shQuote(stats_2), "2>&1"
    ) 
    
    shortstack_cmd_2 <- paste(shortstack_cmd_2,collapse = " ")
    shortstack_cmd_2 <- gsub("^ *| *$", "", shortstack_cmd_2)
    
    # Run ShortStack using system command
    system(shortstack_cmd_2, intern=FALSE)
  }
  
  #####  clean up 
  if (tidy) {
    # List all directories in path_2 and exclude path_2 itself
    map_2_files <- list.dirs(path_2, full.names = TRUE, recursive = TRUE)
    map_2_files <- map_2_files[!map_2_files == path_2]
    
    # Iterate over each directory and remove files that don't match the conditions
    for (k in seq_along(map_2_files)) {
      rm_cmd_2 <- paste0("find '", map_2_files[k], 
     "' -type f ! \\( -name '*.bam' -o -name 'Results.txt' \\) -exec rm {} \\;")
      system(rm_cmd_2, intern=FALSE)
    }
  }
  cat("\n")
  cat("\n")
  message(" --- Mapping of sRNAseq samples is complete --- ")
  message("Results saved to: ", path_1, " & ", path_2)
  message("Loci file saved to: ", loci_out)
}







################## mRNA map  #####################################

mRNA_map <- function(sampleData,
                     input_files_dir, 
                     output_dir,
                     genomefile, 
                     condaenv,
                     annotationfile,
                     threads, 
                     mmap,
                     order ,
                     a, 
                     stranded,
                     mode,
                     nonunique,
                     type,
                     idattr){
  # generate output folders (check if they exist )
  path_1 <- file.path(output_dir, "1_mRNA_preprocessing")
  if (!dir.exists(path_1)) {
    dir.create(path_1)
  }
  
  # check for genome index - genomefile
  stats <- file.path(path_1, "alignment_stats.txt")
  system(paste0(">> ", stats))
  
  # 1 - bowtie - build 
  index_extension <- "ht2"
  base_genomefile <- tools::file_path_sans_ext(basename(genomefile))
  dir_genomfile <- dirname(genomefile)
  files_dir <- c("ls", dir_genomfile)
  files_dir <- paste(files_dir,collapse = " ")
  files_dir <- gsub("^ *| *$", "", files_dir)
  files_dir_res <- system(files_dir, intern=TRUE)
  
  index_found <- grep(paste("^", base_genomefile, ".*\\.", index_extension, 
                            "$", sep = ""), files_dir_res)
  base_genomefile_path <- tools::file_path_sans_ext(genomefile)
  if (length(index_found) > 0) {
    message("Genome index already built.")
  } else {
    message("Genome index not already built. Building with Bowtie... ")
    
    hisat_build_cmd <- c("hisat2-build --threads", threads, 
                          genomefile, base_genomefile_path, ">>", 
                          shQuote(stats), "2>&1") 
    hisat_build_cmd <- paste(hisat_build_cmd,collapse = " ")
    hisat_build_cmd <- gsub("^ *| *$", "", hisat_build_cmd)
    system(hisat_build_cmd, intern=FALSE) ## -- get it to print to summary doc. 
    cat("\n")
    message("Genome index build. ")
  }
  nline <- paste("echo '' >> ",  shQuote(stats))
  system(nline, intern=FALSE)
  
  # paired end 
  if(ncol(sampleData)>2){
    pairend_cmd_message <-paste("echo Running HISAT2 Command: 
    'hist2 -p", shQuote(threads), 
    "-x [genomefile] 
     -1 shQuote([fastqfile_1.fq]) 
     -2 shQuote([fastqfile_2.fq]) 
    | samtools view -bS | samtools sort -o [outputfile.bam]' >>",shQuote(stats))
    
    pairend_cmd_message <- paste(pairend_cmd_message,collapse = " ")
    pairend_cmd_message <- gsub("^ *| *$", "", pairend_cmd_message)
    system(pairend_cmd_message, intern=FALSE)
    nline <- paste("echo '' >> ",  shQuote(stats))
    system(nline, intern=FALSE)
    
    pairend_cmd_2_message <-paste("echo Running HTSeq Command: 
    'python -m HTSeq.scripts.count
    --format=bam 
    --order=", shQuote(order),
    "--a=", shQuote(a),                          
    "--stranded=",shQuote(stranded), 
    "--mode=",shQuote(mode), 
    "--nonunique=",shQuote(nonunique), 
    "--type=",shQuote(type),
    "--idattr=", shQuote(idattr), 
    "[outputfile.bam] [annotationfile.gff]' >>", shQuote(stats))
  
    
    pairend_cmd_2_message <- paste(pairend_cmd_2_message,collapse = " ")
    pairend_cmd_2_message <- gsub("^ *| *$", "", pairend_cmd_2_message)
    system(pairend_cmd_2_message, intern=FALSE)
    nline <- paste("echo '' >> ",  shQuote(stats))
    system(nline, intern=FALSE)
    
      for(i in seq_len(nrow(sampleData))){
        pair_files <- as.character(sampleData[i,])
        sample_name <- pair_files[1]
        fastq_1 <- pair_files[2]
        fastq_2 <- pair_files[3]
        unqiuefolder <- file.path(output_dir,sample_name)
        unqiuefolder_mkdir <- dir.create(unqiuefolder, recursive = TRUE) 
        bam <-file.path(unqiuefolder, paste0(sample_name,".bam"))
        time_cmd <- paste(c("echo", shQuote(as.character(Sys.time())), ">>", 
                            shQuote(stats)), collapse = " ")
        time_cmd <- gsub("^ *| *$", "", time_cmd)
        system(time_cmd, intern=FALSE)
        sample_cmd <- paste(c("echo Working with sample:", 
                              shQuote(sample_name), ">>", shQuote(stats)),
                            collapse = " ")
        sample_cmd <- gsub("^ *| *$", "", sample_cmd)
        system(sample_cmd, intern=FALSE)
        stateline <- paste("echo 1.Alignment: >> ",  shQuote(stats))
        system(stateline, intern=FALSE)

        pairEndmap_cmd <- c("hisat2 -p", threads, 
                        "-x", shQuote(base_genomefile_path),"-1", 
                        shQuote(file.path(input_files_dir, fastq_1)), 
                        "-2",shQuote(file.path(input_files_dir, fastq_2)),
                        "| samtools view -bS | samtools sort -o", 
                        shQuote(bam), "2>>", shQuote(stats))
        pairEndmap_cmd <- paste(pairEndmap_cmd,collapse = " ")
        pairEndmap_cmd <- gsub("^ *| *$", "", pairEndmap_cmd)
        system(pairEndmap_cmd, intern=FALSE)
        stateline <- paste("echo 2.Raw count: >> ",  shQuote(stats))
        system(stateline, intern=FALSE)
        # if unique only:
        if(mmap == "n"){
          out <- file.path(unqiuefolder, paste0(sample_name,"_uniqueReads1.bam"))
          unique_reads <- c("samtools view", shQuote(bam), "| grep", 
                            c("NH:i:1"),">", shQuote(out))
          unique_reads <- paste(unique_reads,collapse = " ")
          unique_reads <- gsub("^ *| *$", "", unique_reads)
          system(unique_reads, intern=FALSE)
          
          # reheader
          header <- file.path(unqiuefolder, paste0(sample_name,"_header.txt"))
          reheader1 <- c("samtools view -H", shQuote(bam), ">", shQuote(header)) 
          reheader1 <- paste(reheader1,collapse = " ")
          reheader1 <- gsub("^ *| *$", "", reheader1)
          system(reheader1, intern=FALSE)
          
          bam_rce <- file.path(unqiuefolder, paste0(sample_name,"_uniqueReads.bam"))
          reheader2 <- c("cat", shQuote(header), shQuote(out), ">",  shQuote(bam_rce))
          reheader2 <- paste(reheader2,collapse = " ")
          reheader2 <- gsub("^ *| *$", "", reheader2)
          system(reheader2, intern=FALSE)
          
          # remove truncated: 
          rmtrunc <- c("rm", shQuote(out))
          rmtrunc <- paste(rmtrunc,collapse = " ")
          rmtrunc <- gsub("^ *| *$", "", rmtrunc)
          system(rmtrunc, intern=FALSE)
          
        } else 
          if (mmap != "n"){
            bam_rce <- bam
          }
        # counts --- 
        count_file <-file.path(unqiuefolder, "Results.txt")
        HTseq_cmd <- c("python -m HTSeq.scripts.count", " ",
                       "--format=bam", " ",
                       "--order=",shQuote(order), " ",
                       "--a=", shQuote(a)," ",
                       "--stranded=",shQuote(stranded)," ",
                       "--mode=",shQuote(mode), " ",
                       "--nonunique=",shQuote(nonunique)," ",
                       "--type=",shQuote(type)," ",
                       "--idattr=", shQuote(idattr), " ",
                       shQuote(bam_rce), " ",
                       shQuote(annotationfile), " ",
                       " > ", shQuote(count_file),
                       " 2>> ", shQuote(stats))
        HTseq_cmd <- paste(HTseq_cmd,collapse = "")
        HTseq_cmd <- gsub("^ *| *$", "", HTseq_cmd)
        system(HTseq_cmd, intern=FALSE)
        
        complete_cmd <- paste(c("echo Complete.", ">>", shQuote(stats)),
                            collapse = " ")
        complete_cmd <- gsub("^ *| *$", "", complete_cmd)
        system(complete_cmd, intern=FALSE)
        nline <- paste("echo '' >> ",  shQuote(stats))
        system(nline, intern=FALSE)
      }
  } else {
    singleend_cmd_message <- paste("echo Running HISAT2 Command: 'hist2 -p", 
    shQuote(threads), "-x [genomefile] 
                    -U shQuote([fastqfile_1.fastq]) | samtools view -bS",
                    "| samtools sort -o [outputfile.bam]' >>", shQuote(stats))
    
    singleend_cmd_message <- paste(singleend_cmd_message,collapse = " ")
    singleend_cmd_message <- gsub("^ *| *$", "", singleend_cmd_message)
    system(singleend_cmd_message, intern=FALSE)
    nline <- paste("echo '' >> ",  shQuote(stats))
    system(nline, intern=FALSE)
    
    singleend_cmd_2_message <- paste(
      "echo 'Running HTSeq Command: python -m HTSeq.scripts.count --format=bam 
      --order=",
      shQuote(order), " ",
      "--a=", shQuote(a), " ",                              
      "--stranded=", shQuote(stranded), " ",
      "--mode=", shQuote(mode), " ",
      "--nonunique=", shQuote(nonunique), " ",
      "--type=", shQuote(type), " ",
      "--idattr=", shQuote(idattr), " ",
      "[outputfile.bam] [annotationfile.gff]' >> ", shQuote(stats)
    )
  
    singleend_cmd_2_message <- paste(singleend_cmd_2_message, collapse = "")
    singleend_cmd_2_message <- gsub("^ *| *$", "", singleend_cmd_2_message)
    system(singleend_cmd_2_message, intern=FALSE)
    
    nline <- paste("echo '' >> ",  shQuote(stats))
    system(nline, intern=FALSE)
    
   for(k in seq_len(nrow(sampleData))) {
     pair_files <- as.character(sampleData[k,])
     sample_name <- pair_files[1]
     fastq_1 <- pair_files[2]
     unqiuefolder <- file.path(path_1,sample_name)
     unqiuefolder_mkdir <- dir.create(unqiuefolder, recursive = TRUE) 
     bam <-file.path(unqiuefolder, paste0(sample_name,".bam"))
     time_cmd <- paste(c("echo", shQuote(as.character(Sys.time())), ">>", 
                         shQuote(stats)), collapse = " ")
     time_cmd <- gsub("^ *| *$", "", time_cmd)
     system(time_cmd, intern=FALSE)
     sample_cmd <- paste(c("echo Working with sample:", 
                           shQuote(sample_name), ">>", shQuote(stats)),
                         collapse = " ")
     sample_cmd <- gsub("^ *| *$", "", sample_cmd)
     system(sample_cmd, intern=FALSE)
     stateline <- paste("echo 1.Alignment: >> ",  shQuote(stats))
     system(stateline, intern=FALSE)
     singleEndmap <- c("(hisat2 -p", threads, 
                     "-x", shQuote(base_genomefile_path),"-U", 
                     shQuote(file.path(input_files_dir,fastq_1)), 
                     "| samtools view -bS | samtools sort -o", 
                     shQuote(bam),") 2>>", shQuote(stats))
     singleEndmap <- paste(singleEndmap,collapse = " ")
     singleEndmap <- gsub("^ *| *$", "", singleEndmap)
     system(singleEndmap, intern=FALSE)
     stateline <- paste("echo 2.Raw count: >> ",  shQuote(stats))
     system(stateline, intern=FALSE)
     
     if(mmap == "n"){
       out <- file.path(unqiuefolder, paste0(sample_name,"_uniqueReads1.bam"))
       unique_reads <- c("samtools view", shQuote(bam), "| grep", 
                         c("NH:i:1"),">", shQuote(out))
       unique_reads <- paste(unique_reads,collapse = " ")
       unique_reads <- gsub("^ *| *$", "", unique_reads)
       system(unique_reads, intern=FALSE)
       
       # reheader
       header <- file.path(unqiuefolder, paste0(sample_name,"_header.txt"))
       reheader1 <- c("samtools view -H", shQuote(bam), ">", shQuote(header)) 
       reheader1 <- paste(reheader1,collapse = " ")
       reheader1 <- gsub("^ *| *$", "", reheader1)
       system(reheader1, intern=FALSE)
       
       bam_rce <- file.path(unqiuefolder, paste0(sample_name,"_uniqueReads.bam"))
       reheader2 <- c("cat", shQuote(header), shQuote(out), ">",  shQuote(bam_rce))
       reheader2 <- paste(reheader2,collapse = " ")
       reheader2 <- gsub("^ *| *$", "", reheader2)
       system(reheader2, intern=FALSE)
       
       # remove truncated: 
       rmtrunc <- c("rm", shQuote(out))
       rmtrunc <- paste(rmtrunc,collapse = " ")
       rmtrunc <- gsub("^ *| *$", "", rmtrunc)
       system(rmtrunc, intern=FALSE)
       
     } else 
       if (mmap != "n"){
         bam_rce <- bam
       }
     # counts --- 
     count_file <-file.path(unqiuefolder, "Results.txt")
     HTseq_cmd <- c("python -m HTSeq.scripts.count"," ", 
                    "--format=bam", " ", 
                    "--stranded=",shQuote(stranded), " ", 
                    "-a=", shQuote(a), " ", 
                    "--mode=",shQuote(mode), " ", 
                    "--nonunique=",shQuote(nonunique), " ", 
                    "--type=",shQuote(type)," ", 
                    "--idattr=", shQuote(idattr)," ", 
                    shQuote(bam_rce), " ",  shQuote(annotationfile)," ",  
                    " > ", shQuote(count_file),  " 2>> ", shQuote(stats))
     HTseq_cmd <- paste(HTseq_cmd,collapse = "")
     HTseq_cmd <- gsub("^ *| *$", "", HTseq_cmd)
     system(HTseq_cmd, intern=FALSE)
     
     complete_cmd <- paste(c("echo Complete.", ">>", shQuote(stats)),
                           collapse = " ")
     complete_cmd <- gsub("^ *| *$", "", complete_cmd)
     system(complete_cmd, intern=FALSE)
     nline <- paste("echo '' >> ",  shQuote(stats))
     system(nline, intern=FALSE)
     }
  }
  cat("\n")
  cat("\n")
  message(" --- Mapping of mRNAseq samples is complete --- ")
  message("Results saved to: ", path_1)
}



################## mRNA map extras  #####################################

exists_conda <- function(package_name){
  # Replace "your_package_name" with the name of the package you want to check
  cmd <- paste(c("conda list | grep", shQuote(package_name)), collapse = " ")
  cmd_b <- gsub("^ *| *$", "", cmd)
  out <- system(cmd_b, intern=TRUE )
  
  if (length(out) == 0) {
    res <- FALSE
  } else {
    res <- TRUE
  }
  return(res)
}

exists_htseq <- function(package_name){
  # Replace "your_package_name" with the name of the package you want to check
  cmd <- paste("htseq-count --version", collapse = " ")
  cmd_b <- gsub("^ *| *$", "", cmd)
  out <- system(cmd_b, intern=TRUE )
  
  if (length(out) == 0) {
    res <- FALSE
  } else {
    res <- TRUE
  }
  return(res)
}



################## global variable storage #####################################
utils::globalVariables(c("ID", "DicerConsensus", "nt_20", "nt_21", "nt_22",
                         "nt_23", "nt_24", "group", "qc", "score", "phase",
                         "type", "V1", "V10", "V11", "V12", "V13", "V14", "V15",
                         "V16", "V17", "V2", "V3", "V4", "V5", "V6", "V7", "V8",
                         "V9", "padj",".", "nt_N", "n", "FPKM_mean", "chr",
                         "key", "Count", "Class", "padjusted", "pvalue", "freq",
                         "value" , "variable" , "repeats_info" , "Genome" ,
                         "Dataset" ,"setNames" , "DicerCall" , "Reads" , "RPM" ,
                         "MajorRNA", "i", "other", "report", "DicerCounts", 
                         "Sequence", "new_df",  "PC1", "PC2", "Conditions",
                         "name", "Length", "Locus", "Name", "SampleCounts", 
                         "UniqueReads", "groups", "is", "log2FoldChange","mRNA",
                         "no", "none", "pos", "significance", "tidy", "a"))

