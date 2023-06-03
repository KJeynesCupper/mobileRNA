#' Import and organise sRNAseq data set
#'
#' @description Organises the sRNAseq results data, outputted by ShorStack, for
#' each sample replicate into a single dataframe using the loci information
#' produced by the [RNAlocate::RNAloci()] function which contains a
#' list of all the clusters identified across samples.
#'
#'
#' @details
#' Supply the directory to the folder storing the various sample replicate folders
#' produced by ShortStack analysis. Within each sample replicate folder, there is
#' "Result.txt". This file contains many different columns, but
#' for the analysis the columns of interest are  "Locus", "DicerCall", "Reads", and
#' "RPM" which will all be imported and loads into the new data frame. Note that
#' "Reads" will be converted to "Counts" as output of the function.
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
#' Counts, MajorRNA and RPM. The DicerCall represents the size of most
#' abundant small RNA size based on the parameter used in ShortStack.
#' The Count column represents the number of aligned sRNA-seq reads that overlap
#' the locus. The RPM represents the reads per million and rhe MajorRNA represents
#' the most common RNA sequence.
#' For each replicate included in the analysis, these columns are labeled with
#' the type and then then name of the sample, for example, for
#' a sample called "Sample1", the information from this sample will be stored in
#' columns DicerCall_Sample1, Count_Sample1 and RPM_Sample1.
#'
#'
#'
#' @param loci dataframe; rows represent individual sRNA dicer-derived clusters
#' columns state genomic location (labelled "Locus") and the cluster name
#' (labeled "Cluster"). This file contains an accumulation of all
#' sRNA dicer-derived clusters  identified across samples in the analyse, without
#' duplication. This can be produced using the [RNAlocate::RNAloci()]
#' function
#'
#'
#'
#' @param directory Path to directory containing of sample folders. NOTE: Following
#' the suggested pre-processing steps, these can be found in second alignment
#' folders.
#'
#'
#' @param samples Vector of characters naming the sample names correlating
#' to outputted folders located in the \code{`directory`} argument path.
#'
#'
#' @param report Logical; prompts progress.
#'
#'
#'
#'
#'
#' @examples
#' \dontrun{
#'
#' clusters <- utils::read.table(file = "./data/reference/ClustersInfo.txt",
#' header = TRUE, sep = "\t", stringsAsFactors = TRUE,comment.char="")
#'
#'RNAorganise(loci = clusters,
#'            directory = "./analysis/alignment_unique_two/",
#'            samples = c( "TomEgg_1", "TomEgg_2", "TomEgg_3",
#'                        "TomTom_1" , "TomTom_2" , "TomTom_3"))
#'
#'
#'}
#'
#' @export
#' @importFrom data.table "data.table"
#' @importFrom data.table "setnames"
#' @importFrom  utils "read.table"
#' @importFrom data.table "fread"
#' @importFrom dplyr "mutate"
#' @importFrom dplyr "across"
#' @importFrom dplyr "contains"
#' @importFrom tidyr "replace_na"
#' @importFrom data.table ":="
RNAimport <- function(loci, directory, samples,
                        report = TRUE) {

  if (base::missing(loci)) {
    stop("loci must be an object of data frame or directory to plain text dataframe.
         Dataframe represents the loci of clusters identified across all sample replicates")
  }

  # import loci data
  if(base::inherits(loci, c("data.frame")) == FALSE){
    message("Loading loci file from directory...")
    loci <- utils::read.table(file = loci,header = TRUE, sep = "\t",
                              stringsAsFactors = TRUE,comment.char="")
  }
  # check loci file has column called locus
  col_needed_loci <- "Locus"

  if (length(setdiff(col_needed_loci, names(loci))) > 0 || grep("Locus", colnames(loci)) != 1 ) {
    stop("loci dataframe must contain 'Locus' column,  each row contains the
         genomic coordinates for a given cluster.

         'Locus' column must  be in the first position within the dataframe ")
  }

  # LOad sample data as list of data frames, with index as file name.
  dt_list <- list()
  for (file in samples) {
    dt_list[[file]] <- data.table::fread(paste0(directory, file, "/Results.txt"), header = TRUE)
  }
  # remove any hashtags from header (shortstack add this to header line, position 1)
  dt_list <- lapply(dt_list, function(x) setNames(x, gsub("#", "", names(x))))
  # Check each data frame in the list for the required columns
  required_cols <- c("Locus", "DicerCall", "Reads", "RPM", "MajorRNA")
  for (df in dt_list) {
    if (!all(required_cols %in% colnames(df))) {
      stop("Sample data frame does not contain all required columns: ",
           paste(setdiff(required_cols, colnames(df)), ".", collapse = ", ",
                 "Make sure there is not a hashtag or similar in the header line
                 of the input file(s)"))
    }
  }


  # Convert all input dataframes to data.tables
  dt1 <- data.table::data.table(loci)

  # create list of data to add, name each with its proper name.
  #dt_list <- list(...)
  #names(dt_list) <- as.character(substitute(list(...))[-1])
  #dt_list <- lapply(dt_list, data.table::data.table)

  # Define a function to update df1 with the matching values from a single input dataframe
  update_locus_df <- function(dt, i) {
    # Join df1 and the current input dataframe on chromosome and coordinate range
    join_cols <- c("Locus")
    dt_match <- dt1[dt, on = join_cols, nomatch = 0]

    # Aggregate the matching rows by chromosome, start coordinate, and end coordinate,
    # and compute the sum of DicerCall, Reads, and RPM values for each group
    dt_agg <- dt_match[, .(DicerCall = as.character(DicerCall),
                           Count=sum(Reads),
                           RPM = sum(RPM)),
                       MajorRNA = MajorRNA,
                       by = join_cols]

    # Rename the aggregated columns
    col_names <- paste0(c("DicerCall_", "Count_", "RPM_", "MajorRNA_"),  i)
    data.table::setnames(dt_agg, c("Locus", col_names))

    # Merge the aggregated values back into df1
    dt1[dt_agg, on = join_cols, (col_names) := mget(col_names)]

    # Print progress
    if (report) {
      cat("Added information from", i, "to the analysis locus dataframe" ,"\n")
    }
  }

  # Update df1 with the matching values from each input dataframe
  for (i in seq_along(dt_list)) {
    update_locus_df(dt_list[[i]], names(dt_list)[i])
  }

  # Fill in missing values with 0 or N
  ## Dicer call needs to character/factor
  dt1 <- dt1 %>%
    dplyr::mutate(dplyr::across(dplyr::contains('Count_'), ~tidyr::replace_na(.,0))) %>%
    dplyr::mutate(dplyr::across(dplyr::contains('RPM_'), ~tidyr::replace_na(.,0))) %>%
    dplyr::mutate(dplyr::across(dplyr::contains('MajorRNA_'), ~tidyr::replace_na(.,0))) %>%
    dplyr::mutate(dplyr::across(dplyr::contains('DicerCall_'), ~tidyr::replace_na(.,"N")))


  # Convert dt1 back to a data.frame and return it
  res_data <- as.data.frame(dt1)

  # Split the Locus column into three new columns
  locus_cols <- data.frame(
    chr = sapply(strsplit(res_data$Locus, split = ":"), "[[", 1),
    start = sapply(strsplit(sapply(strsplit(res_data$Locus, split = ":"), "[[", 2), split = "-"), "[[", 1),
    end = sapply(strsplit(sapply(strsplit(res_data$Locus, split = ":"), "[[", 2), split = "-"), "[[", 2)
  )
  df_final <- cbind(res_data[,1], locus_cols, res_data[, 2:ncol(res_data)])
  names(df_final)[1] <- "Locus"


    return(df_final)
}

