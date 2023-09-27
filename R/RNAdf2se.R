#' Convert `mobileRNA` dataframe to SummarizedExperiment object 
#'
#'
#'@description 
#'Convert any `mobileRNA` output dataframe into a SummarizedExperiment object.
#'
#'@param data data.frame produced by the \pkg{mobileRNA} package. 
#'@param input character; must be either "sRNA" or "mRNA"
#'@details 
#'The function relies on the naming structure of columns created by functions
#'in the \pkg{mobileRNA} package, specifically [mobileRNA::RNAimport()] and 
#'[mobileRNA::RNAdicercall()]. 
#'It is able to extract the sample names based on these additions, and organise 
#'the data appropriately. 
#'
#'**For sRNAseq data** 
#'* The rownames contain the locus name and the cluster name. 
#'* The assays represent the additional information including DicerCall, Count, RPM, MajorRNA. 
#'* The rowData includes the Cluster ID, the DicerCounts & the DicerConsensus
#'* The colnames represents the sample replicate names
#'
#'**For mRNAseq data** 
#'* The rownames contain the gene names
#'* The assays represent the additional information including Count &FPKM. 
#'* The rowData includes the gene & the SampleCounts. 
#'* The colnames represents the sample replicate names
#'
#' @return A `SummarizedExperiment` object containing information from working 
#' data frame.  
#'
#'@examples
#' # load data.frame
#'data("sRNA_data")
#'
#'se <- RNAdf2se(input = "sRNA", data = sRNA_data)
#'
#'@export
#'@importFrom SummarizedExperiment "SummarizedExperiment"
#'@importFrom GenomicRanges "GRanges"
#'@importFrom IRanges "IRanges"
#'@importFrom S4Vectors "DataFrame"
#'@importFrom dplyr "select"
RNAdf2se <- function(input= c("sRNA", "mRNA"), data){
  if (base::missing(data)) {
    stop("data is missing. data must be an object of class matrix, data.frame, 
         DataFrame")
  }
  if (!base::inherits(data, c("matrix", "data.frame", "DataFrame"))) {
    stop("data must be an object of class matrix, data.frame, DataFrame")
  }
  if (missing(input) || !is.character(input)) {
    stop("style parameter is missing or not a character vector.")
  }
  
  # Check if input is one of the allowed values
  allowed_input <- c("sRNA", "mRNA")
  if (!input %in% allowed_input) {
    stop("input parameter must be one of 'sRNA' or 'mRNA'.")
  }
  
  if(input == "sRNA"){
    # create gene locations
    if("DicerCounts" %in% colnames(data)) {
      rowRanges <- GenomicRanges::GRanges(data$chr,
                                          IRanges::IRanges(data$start,data$end),
                                          Cluster = data$Cluster, 
                                          DicerCounts = data$DicerCounts,
                                          DicerConsensus= data$DicerConsensus)
    } else {
      rowRanges <- GenomicRanges::GRanges(data$chr,
                                          IRanges::IRanges(data$start,data$end),
                                          Cluster = data$Cluster)
    }
    
    # create table for the columns
    col_names <- colnames(data) # extract sample names 
    sample_names <- unique(sapply(base::strsplit(col_names, "_"), 
                                  function(x) paste(x[-1], collapse = "_"))) 
    sample_names <- sample_names[nzchar(sample_names)] # remove empties
    
    colData <- S4Vectors::DataFrame(samples=sample_names,
                                    row.names=sample_names)
    
    
    # make matrix for all information 
    extra <- S4Vectors::DataFrame(data[, !names(data) %in% c("chr", "start",
                                                             "end", 
                                                             "Cluster", "Locus", 
                                                             "DicerCounts", 
                                                             "DicerConsensus"), 
                                       drop = FALSE])
    extra_vrs <-  gsub("^(.*?)_.*", "\\1", colnames(extra))
    extra_vrs_unique <- unique(extra_vrs)
    element_counts <- table(extra_vrs) # Count occurrences of each element
    occurrences_check <- all(element_counts == length(sample_names)) 
    # Check if all elements have same occurrences
    if(occurrences_check == FALSE){
      stop(paste("data is missing information for all replicates"))
    }
    assay_list <- list()
    # create matrix for each extra data and count  
    for(i in seq_along(extra_vrs_unique) ){
      extra_matrix <-as.matrix(data %>% dplyr::select(
        starts_with(extra_vrs_unique[i])))
      colnames(extra_matrix) <- NULL
      assay_list[[i]] <- extra_matrix
    }
    names(assay_list) <- extra_vrs_unique
    # create SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(assays= assay_list,
                                                     rowRanges=rowRanges, 
                                                     colData=colData)
    rownames(se) <- paste0(data$Locus, " (",data$Cluster, ")" )
  } else 
    if( input == "mRNA"){
      # create gene locations
      
      rowRanges <- GenomicRanges::GRanges(data$chr,
                                          IRanges::IRanges(data$start, 
                                                           data$end, 
                                                           data$width),
                                          Gene = data$Gene, 
                                          SampleCounts = data$SampleCounts)
      
      # create table for the columns
      col_names <- colnames(data) # extract sample names 
      sample_names <- unique(sapply(base::strsplit(col_names, "_"), 
                                    function(x) paste(x[-1], collapse = "_"))) 
      sample_names <- sample_names[nzchar(sample_names)] # remove empties
      
      colData <- S4Vectors::DataFrame(samples=sample_names,
                                      row.names=sample_names)
      
      
      # make matrix for all information 
      extra <- S4Vectors::DataFrame(data[, !names(data) %in% c("chr", "start", 
                                                               "end", 
                                                               "Gene", "Locus",
                                                               "width", 
                                                               "SampleCounts"), 
                                         drop = FALSE])
      extra_vrs <-  gsub("^(.*?)_.*", "\\1", colnames(extra))
      extra_vrs_unique <- unique(extra_vrs)
      element_counts <- table(extra_vrs) # Count occurrences of each element
      occurrences_check <- all(element_counts == length(sample_names)) 
      # Check if all elements have same occurrences
      if(occurrences_check == FALSE){
        stop(paste("data is missing information for all replicates"))
      }
      assay_list <- list()
      # create matrix for each extra data and count  
      for(i in seq_len(nrow(extra_vrs_unique))){
        extra_matrix <-as.matrix(data %>% dplyr::select(
          starts_with(extra_vrs_unique[i])))
        colnames(extra_matrix) <- NULL
        assay_list[[i]] <- extra_matrix
      }
      names(assay_list) <- extra_vrs_unique
      # create SummarizedExperiment object
      se <- SummarizedExperiment::SummarizedExperiment(assays= assay_list,
                                                       rowRanges=rowRanges, 
                                                       colData=colData)
      rownames(se) <- data$Gene 
    }
  return(se)
}
