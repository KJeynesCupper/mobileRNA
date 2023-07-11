#' Import and organise sRNAseq or mRNA data sets
#'
#' @description Load and organise either sRNAseq or mRNAseq results into a
#' single dataframe containing all experimental replicates specified where rows
#' represent either a sRNA locus or gene, respectively.
#'
#'
#' @details
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
#'
#' @param report Logical; prompts progress.
#'
#'@param tidy Logical; removes genes from analysis where there are zero counts
#'across all replicates.
#'
#'
#'@return A dataframe where rows represent sRNA clusters and columns represent
#'replicate information. Replicate information includes Dicercall, Counts, RPM
#'and MajorRNA sequence. Each replicate information is distinguishable as
#'the replicate name is ajoined as a prefix to each column name.
#'
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

RNAimport <- function(input = c("sRNA", "mRNA"), directory, samples,
                        report = TRUE,
                        tidy = TRUE) {

  if(input=="sRNA"){
    # LOad sample data as list of data frames, with index as file name.
    dt_list <- list()
    total_files <- length(samples)
    for (file in samples) {
      options(datatable.showProgress = FALSE)
      dt_list[[file]] <- data.table::fread(paste0(directory, file,
                                                  "/Results.txt"),header = TRUE)
      progress_counter <- file
      progress_message <- paste0("Processing file ", progress_counter, " of ", total_files)
      cat(sprintf("\r%s", progress_message))
      flush.console()
    }
    cat("\n")  # Print a newline after progress is complete
    message("Completed importation of data from directory.")
    
    
    # remove any hashtags from header - added by shortstack
    dt_list <- lapply(dt_list, function(x) setNames(x, gsub("#", "", names(x))))
    # Check each data frame in the list for the required columns
    message("Checking data content...")
    required_cols <- c("Locus", "DicerCall", "Reads", "RPM", "MajorRNA")
    for (df in dt_list) {
      if (!all(required_cols %in% colnames(df))) {
        stop("Sample data frame does not contain all required columns: ",
             paste(setdiff(required_cols, colnames(df)), ".", collapse = ", ",
                   "Make sure there is not a hashtag or similar in the header
                   line of the input file(s)"))
      }
    }
    message("Data content is correct.")
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

      # Print progress
      if (report) {
        cat("Added information from", i, "to the analysis locus dataframe","\n")
      }
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



    # return values
    return(df_final)

  } else
    if(input == "mRNA"){
      # load data as list
      sample_data <- list()
      for (file in samples) {
        sample_data[[file]] <- data.table::fread(paste0(directory, file,
                                                        ".txt"), header = FALSE)
        colnames(sample_data[[file]])[1] <- "Gene"
        colnames(sample_data[[file]])[2] <- "Count"
      }

      # check each file has two columns
      for (df in sample_data) {
        if (!all(ncol(df) == 2)) {
          stop("mRNA dataset(s) does not contain required columns")
        }
      }
      # merge first columns to create list of genes across all samples
      genes <- lapply(sample_data, "[", , "Gene")
      genes_all <- unique(Reduce(merge,genes))

      # ADDs sample information to the genes_all object
      for (i in seq_along(sample_data)){
        matches <- genes_all[sample_data[[i]], on = "Gene", nomatch = 0]
        matches_values <- matches[, .(Count=sum(Count)),by = "Gene"]
        # Rename the aggregated columns
        col_name <- paste0("Count_", names(sample_data)[i])
        data.table::setnames(matches_values, c("Gene", col_name))
        # Merge the aggregated values back into df1
        genes_all[matches_values, on = "Gene", (col_name) := mget(col_name)]
        # Print progress
        if (report) {
          message(paste0("Added information from", i,
                         "to the analysis mRNA dataframe" ,"\n"))
        }
      }

      # Fill in missing values with 0 or N
      mRNA_information <- genes_all %>%
        dplyr::mutate(dplyr::across(dplyr::contains('Count_'),
                                    ~tidyr::replace_na(.,0)))
      # Convert data.frame and return it
      mRNA_information <- as.data.frame(mRNA_information)

      # remove rows with
      if (tidy){
        mRNA_information <- mRNA_information %>%
          dplyr::filter(dplyr::if_any(dplyr::where(is.numeric), ~. != 0))
      }
      return(mRNA_information)
    }
}

