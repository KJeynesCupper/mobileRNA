#' Import and organise sRNAseq or mRNA data sets
#'
#' @description Load and organise either sRNAseq or mRNAseq results into a
#' single dataframe containing all experimental replicates specified where rows
#' represent either a sRNA locus or gene, respectively.
#'
#'
#' @details
#' **For sRNAseq:**
#' Supply the directory to the folder storing the various sample replicate
#' folders produced by ShortStack analysis. Within each sample replicate folder,
#' there is "Result.txt". This file contains many different columns, but
#' for the analysis the columns of interest are  "Locus", "DicerCall", "Reads",
#' and "RPM" which will all be imported and loads into the new data frame.
#' Note that "Reads" will be converted to "Counts" as output of the function.
#' Locus contains the genomic locations of clusters,
#' Dicercall contains the most likely size of the cluster, Reads contains the
#' number of reads mapped to the cluster (recommended to use unique mapping,
#' hence these will be the number of uniquely mapped reads, ie, not including
#' multimapped reads). Lastly, RPM contains the Reads per Million score.
#'
#' \code{"RNAimport"} takes the results from  pre-processed sRNAseq data which
#' has been mapped and undertaken cluster analysis through ShortStack, organises
#'  into a single data frame. Each row represent an sRNA dicer-derived cluster
#'  within the analysis and columns 1 - 5 supply information on the cluster
#'  including the locus, the separated coordinates
#' (chromosome, start, end) and cluster name.
#'
#' Further columns represent data imported for each samples including DicerCall,
#' Counts and RPM. The DicerCall represents the size of most
#' abundant small RNA size based on the parameter used in ShortStack.
#' The Count column represents the number of aligned sRNA-seq reads that overlap
#' the locus. The RPM represents the reads per million.
#' For each replicate included in the analysis, these columns are labeled with
#' the type and then then name of the sample, for example, for
#' a sample called "Sample1", the information from this sample will be stored in
#' columns DicerCall_Sample1, Count_Sample1 and RPM_Sample1.
#'
#'
#'**For mRNAseq:**
#'
#' Supply the directory to the folder storing the various sample replicate
#' files outputted by the alignment program. Each file contains many the gene
#' name (column 1) and the raw count value (column 2). The function identifies 
#' all genes across the analysis, and organises the count values for each sample
#' replicate into a column. For example, this will look like 
#' "Counts_heterograft_1" for the sample replicate with the file name
#' "heterograft_1". Further columns include the gene locus information (chr, 
#' start, end, width) are add to the start of the dataframe, and the FPKM 
#' for each sample replicate is calculated. The FPKM columns are organised and 
#' titled in the same format at the "Counts_" columns. 
#' 
#' When working with mRNAseq data, a GFF genome annotation file is required. 
#' When working with a chimeric biological system and a merged annotation file, 
#' ensure the genomes remain distinguishable. The function 
#' `RNAmergeAnnotations()` can help with this.
#' 
#'@param input string; define type of Next-Generation Sequencing dataset
#'originates from, either "sRNA" or "mRNA" are the only valid inputs.
#'
#' @param directory Path to directory containing of sample folders. NOTE:
#' Following the suggested pre-processing steps, these can be found in second
#' alignment folders.
#'
#'
#' @param samples Vector of characters naming the sample names correlating
#' to outputted folders located in the \code{`directory`} argument path.
#'
#' @param annotation Path to GFF genome annotation file. Required only for 
#' mRNAseq data. 
#'
#'@param tidy Logical; removes genes from analysis where there are zero counts
#'across all replicates.
#'
#'
#'@return 
#'**For sRNAseq:**
#'A dataframe where rows represent sRNA clusters and columns represent
#'replicate information. Replicate information includes Dicercall, Counts, RPM
#'and MajorRNA sequence. Each replicate information is distinguishable as
#'the replicate name is ajoined as a prefix to each column name.
#'
#'**For mRNAseq:**
#'A dataframe where rows represent gene names and columns representing
#'replicate information. Replicate information includes Counts and FPKM. 
#'Additional columns include SampleCounts (the number of replicates which 
#'contains counts), chr (chromosome), start and end coordinates are 
#'also included in the output. Each replicate information is distinguishable as 
#'the replicate name is adjoined as a prefix to each column name.
#' @examples
#' \dontrun{
#'
#' # import sRNAseq data
#' df_sRNA <- RNAimport(input = "sRNA",
#'                      directory = "./analysis/sRNA_mapping_results/",
#'                      samples = c("heterograft_1", "heterograft_2",
#'                      "heterograft_3","selfgraft_1" , "selfgraft_2" ,
#'                      "selfgraft_3"))
#'
#'# import mRNAseq data
#' df_mRNA <- RNAimport(input = "mRNA",
#'                      directory = "./analysis/mRNA_counts/",
#'                      samples = c("heterograft_1", "heterograft_2",
#'                      "heterograft_3","selfgraft_1" , "selfgraft_2" ,
#'                      "selfgraft_3"))
#'
#'}
#'# The output of this function can be explore in the data object sRNA_data
#' data("sRNA_data")
#' head(sRNA_data)
#'
#'
#' @export
#' @importFrom data.table "data.table"
#' @importFrom data.table "setnames"
#' @importFrom utils "read.table"
#' @importFrom data.table "fread"
#' @importFrom dplyr "mutate"
#' @importFrom dplyr "across"
#' @importFrom dplyr "contains"
#' @importFrom tidyr "replace_na"
#' @importFrom data.table ":="
#' @importFrom magrittr "%>%"
#' @importFrom dplyr "filter"
#' @importFrom dplyr "if_any"
#' @importFrom dplyr "where"
#' @importFrom stats "setNames"
#' @importFrom utils "flush.console"
#' @importFrom stats "complete.cases"
#' @importFrom rtracklayer "import.gff"
#' @importFrom S4Vectors "mcols"
#' @importFrom Repitools "annoGR2DF"
#' @importFrom dplyr "select"
#' @importFrom dplyr "rename"
#' @importFrom dplyr "rename_with"
#' @importFrom stringr "str_detect"
RNAimport <- function(input = c("sRNA", "mRNA"), directory, samples,
                      annotation = NULL, tidy = TRUE) {
  if (base::missing(input) || !input %in% c("sRNA", "mRNA")) {
    stop(paste("Please the data-type to the `input` paramter"))
  }
  if (base::missing(directory)) {
    stop(paste("Please specify a accessable directory where files are stored"))
  }
  if (base::missing(samples)) {
    stop(paste("Please specify a vector storing sample names matching files"))
  }
  
  if(input=="sRNA"){
    # LOad sample data as list of data frames, with index as file name.
    dt_list <- list()
    total_files <- length(samples)
    file_n <- 0
    for (file in samples) {
      file_n <- file_n + 1
      options(datatable.showProgress = FALSE)
      dt_list[[file]] <- data.table::fread(paste0(directory, file,
                                                  "/Results.txt"),header = TRUE)
      progress_counter <- file
      progress_message <- paste0("Processing sample: ", progress_counter,".", "\n",
                                 "File ", file_n, " of ", total_files)
      cat(sprintf("\r%s", progress_message))
      utils::flush.console()
    }
    cat("\n")  # Print a newline after progress is complete
    cat("Completed importation of data from directory. \n")
    
    
    # remove any hashtags from header - added by shortstack
    dt_list <- lapply(dt_list, function(x) setNames(x, gsub("#", "", names(x))))
    # Check each data frame in the list for the required columns
    cat("Checking data content... \n")
    required_cols <- c("Locus", "DicerCall", "Reads", "RPM", "MajorRNA")
    for (df in dt_list) {
      if (!all(required_cols %in% colnames(df))) {
        stop("Sample data frame does not contain all required columns: ",
             paste(setdiff(required_cols, colnames(df)), ".", collapse = ", ",
                   "Make sure there is not a hashtag or similar in the header
                   line of the input file(s)"))
      }
    }
    cat("Data content is correct. \n")
    cat("\n") 
    
    # merge first columns to create list of loci across all samples
    loci <- lapply(dt_list, "[", , "Locus")
    loci_all <- unique(Reduce(merge,loci))

    # Define a function to update the loci with the matching values from a
    # single input dataframe
    update_locus_df <- function(dt, i) {
      # Join loci and the current input df on chromosome and coordinate range
      join_cols <- c("Locus")
      dt_match <- loci_all[dt, on = join_cols, nomatch = 0]

      # Aggregate the matching rows by chromosome, start coordinate, & end coor,
      # and compute the sum of DicerCall, Reads, and RPM values for each group
      dt_agg <- dt_match[, .(DicerCall = as.character(DicerCall),
                             Count=sum(Reads),
                             RPM = sum(RPM),
                             MajorRNA = as.character(MajorRNA)),
                         by = join_cols]

      # Rename the aggregated columns
      col_names <- paste0(c("DicerCall_", "Count_", "RPM_", "MajorRNA_"),  i)
      data.table::setnames(dt_agg, c("Locus", col_names))

      # Merge the aggregated values back into df1
      loci_all[dt_agg, on = join_cols, (col_names) := mget(col_names)]

    }

    # Update loci with the matching values from each input dataframe
    for (i in seq_along(dt_list)) {
      update_locus_df(dt_list[[i]], names(dt_list)[i])
    }

    # Fill in missing values with 0 or N
    ## Dicer call needs to character/factor
    loci_all <- loci_all %>%
      dplyr::mutate(dplyr::across(dplyr::contains('Count_'),
                                  ~tidyr::replace_na(.,0))) %>%
      dplyr::mutate(dplyr::across(dplyr::contains('RPM_'),
                                  ~tidyr::replace_na(.,0))) %>%
      dplyr::mutate(dplyr::across(dplyr::contains('DicerCall_'),
                                  ~tidyr::replace_na(.,"N")))%>%
      dplyr::mutate(dplyr::across(dplyr::contains('MajorRNA_'),
                                  ~tidyr::replace_na(.,"N"))) %>%
      dplyr::mutate_all(~ ifelse(. == "*", "N", .)) # remove any astriks to "N"


    # Convert loci_all back to a data.frame and return it
    res_data <- as.data.frame(loci_all)

    # Split the Locus column into three new columns
    locus_cols <- data.frame(
      chr = sapply(strsplit(res_data$Locus, split = ":"), "[[", 1),
      start = sapply(strsplit(sapply(strsplit(res_data$Locus, split = ":"),
                                     "[[", 2), split = "-"), "[[", 1),
      end = sapply(strsplit(sapply(strsplit(res_data$Locus, split = ":"),
                                   "[[", 2), split = "-"), "[[", 2)
    )
    df_final <- cbind(res_data[,1], locus_cols, res_data[, 2:ncol(res_data)])
    names(df_final)[1] <- "Locus"
    # order by chr
    df_final <- df_final[order(df_final$chr),]
    # insert cluster name
    cluster_names <-  paste0("cluster_", 1:nrow(df_final))
    df_final <- as.data.frame(append(df_final, list(Cluster = cluster_names),
                                     after = 4))
    
    # Remove rows with no counts 
    count_columns <- grep("^Count", names(df_final))
    # Identify rows where all values in Count columns are zero
    rows_to_remove <- apply(df_final[count_columns], 1, function(row) all(row == 0))
    # Remove rows with all zero values in Count columns
    df_final <- df_final[!rows_to_remove, ]   
    # return values
    return(df_final)

  } else
    if(input == "mRNA"){
      if (base::missing(annotation)) {
        stop(paste("Please specify a accessable path to a GFF file"))
      }
      # load data as list
      sample_data <- list()
      for (file in samples) {
        sample_data[[file]] <- data.table::fread(paste0(directory, file,
                                                        ".txt"), header = FALSE)
        colnames(sample_data[[file]])[1] <- "Gene"
        colnames(sample_data[[file]])[2] <- "Count"
     
      }
      # remove rows with extra info 
      sample_data <- lapply(sample_data, function(x) {
        x[!grepl("__", x$Gene),]
      })
      
      
      # check each file has two columns
      for (df in sample_data) {
        if (!all(ncol(df) == 2)) {
          stop("mRNA dataset(s) does not contain required columns")
        }
      }
      # merge first columns to create list of genes across all samples
      genes <- lapply(sample_data, "[", , "Gene")
      genes_all <- unique(Reduce(merge,genes))
      
      # add gene length
      annotation_file <- rtracklayer::import.gff(annotation)
      genes_info <- annotation_file[which(S4Vectors::mcols(annotation_file)$type == "gene")]
      genes_info <- Repitools::annoGR2DF(genes_info)  %>% 
        dplyr::select(chr, start, end, Name, width)%>% 
        dplyr::rename(Gene = Name)
      Locus <- paste0(genes_info$chr, ":",genes_info$start,"-",
                            genes_info$end)
      genes_info <- cbind(Locus, genes_info)
  
      # merge gene list with annotation info. 
      merged_gene_info <- merge(genes_all, genes_info, by = "Gene", all.x = TRUE)
      #merged_gene_info <- merged_gene_info[stats::complete.cases(merged_gene_info), ]
      gene_widths <- merged_gene_info$width
      
      # ADDs sample information to the genes_all object
      for (i in seq_along(sample_data)){
        matches <- merged_gene_info[sample_data[[i]], on = "Gene", nomatch = 0]
        matches_values <- matches[, .(Count=sum(Count)),by = "Gene"]
        
        # Rename the aggregated columns
        col_name <- paste0("Count_", names(sample_data)[i])
        data.table::setnames(matches_values, c("Gene", col_name))
        # Merge the aggregated values back into df1
        merged_gene_info[matches_values, on = "Gene", (col_name) := mget(col_name)]
        # Print progress
        if (report) {
          cat(paste0("Added information from ", names(sample_data[i]) ," to the analysis mRNA dataframe" ,"\n"))
        }
      }
      
      # Fill in missing values with 0 or N
      mRNA_information <- merged_gene_info %>%
        dplyr::mutate(dplyr::across(dplyr::contains('Count_'),
                                    ~tidyr::replace_na(.,0)))
      # set genes as rownames 
      fpkm <- apply(X = subset(mRNA_information, 
                               select = c(-Gene, -chr, -start, -end, -width)),
                    MARGIN = 2, 
                    FUN = function(x) {
                      sum_x <- sum(as.numeric(x))
                      if (sum_x == 0) {
                        t <- 0
                      } else {
                        t <- 10^9 * x / gene_widths / sum_x
                      }
                      t
                    })
      # add prefix to all columns 
     fpkm <- data.frame(fpkm)
     names_col <- sub("^Count_", "", colnames(fpkm))
     colnames(fpkm) <- paste0("FPKM_", names_col)
     
      
      # add result to df to output
      mRNA_information <- cbind(mRNA_information, fpkm)
      
      # add thresholding values - number of replicates with counts 
      t <- names(mRNA_information)[grep("^Count_", names(mRNA_information))]
      SampleCounts_vals <- apply(X = subset(mRNA_information, 
                               select = c(t)),
                    MARGIN = 1, 
                    FUN = function(x) {
                      length(x[x > 0])
                    })
      
      mRNA_information$SampleCounts <- SampleCounts_val
      
      # remove rows with
      if (tidy){
        mRNA_information <- mRNA_information %>%
          dplyr::filter(SampleCounts != 0)
      }
      return(mRNA_information)
    }
}

