#' Import and organise sRNAseq data set
#'
#' @description Imports sRNAseq pre-processed data outputted by ShorStack and
#' organises into a single data frame.
#'
#' @details
#' \code{"RNAimport"} takes the results from  pre-processed sRNAseq data which
#' has been mapped and undertaken cluster analysis through ShortStack, and
#' imports it from a given directory. Data is organised into a single data frame
#' where each row represent an sRNA dicer-derived cluster. Columns 1 - 5
#' supply information on the cluster including the name, coordinates
#' (chromosome, start, end) and width of locus.
#'
#' Further columns represent data imported for each samples including cluster
#' name, DicerCall, Counts and RPM. The DicerCall represents the size of most
#' abundant small RNA size based on the parameter used in ShortStack.
#' The Count column represents the number of aligned sRNA-seq reads that overlap
#' the locus. The RPM represents the reads per million.
#' For each replicate included in the analysis, these columns are labeled with
#' the type and then then name of the samples, for example, DicerCall_Sample1.
#'
#'
#' To import the samples, the `samples` argument must represent a vector of
#' the names of the folders present in the directory containing the sample
#' replicate, ie, replicate names. It is recommended to use results from the
#' unique alignment found in \code{"./second_alignment_unique/"} over multimapped
#' reads. It is recommended to run a parallel analysis one using the uniquely
#' mapped reads and the other using multimapped reads.
#'
#' @param results  Path to directory containing of sample folders. NOTE: Following
#' the suggested pre-processing steps, these can be found in second alignment
#' folders.
#'
#' @param samples Vector of characters naming the sample names correlating
#' to outputted folders either \code{"second_alignment_unique"} or
#' \code{"second_alignment_multi"}.
#'
#' @param clusters Import \code{.gff3} from \code{"reference"} folder.
#'
#'
#' @param features Path to a .bed file containing the genomic information
#' (annotations/features) which overlaps with the locations of the
#' dicer-derived clusters. The file is created as a part of the pre-processing
#' steps, see [RNAlocate::findOverlap()] function for more information. This is
#' an optional argument, but will streamline your analysis later on if you wish
#' to identify the genomic location of the cluster and the type of feature it
#' overlaps with.
#'
#'
#'
#' @examples
#' \dontrun{
#' genomic_features <- "./data/reference/overlap.bed"
#' directory <- "./data/alignment_unique_two/"
#' sample_names <- c("TomEgg_1","TomEgg_2", "TomEgg_3",
#'               "TomTom_1", "TomTom_2", "TomTom_3" )
#'
#' clusters <- rtracklayer::import.gff("./data/reference/ClustersInfo.gff3")
#'
#' sRNA_data <- RNAimport(results = directory,
#'                            samples = sample_names,
#'                            clusters = clusters,
#'                            features = genomic_features)
#'}
#'
#' @export
#' @importFrom utils "read.table"
#' @importFrom Repitools "annoGR2DF"
#' @importFrom dplyr "%>%"
#' @importFrom S4Vectors "mcols"
#' @importFrom BiocGenerics "width"
#' @importFrom dplyr "select"
#' @importFrom dplyr "left_join"
RNAimport <- function(results, samples, clusters, features = NULL){
  if (base::missing(results)|| !base::inherits(results, c("character"))) {
    stop("results must be an object of character, representing a path to a
         directory")
  }
  if (base::missing(samples)|| !base::inherits(samples, c("character"))) {
    stop("samples must be an object of character, representing a vector of
         sample names")
  }
  if (base::missing(clusters)|| !base::inherits(clusters, c("GRanges"))) {
    stop("clusters must be an object of GRanges")
  }

  for (i in samples) {
    array <- utils::read.table(paste0(results,i,"/Results.txt", sep=""),
                               sep="\t", header =F, row.names=NULL)
    suppressWarnings(S4Vectors::mcols(clusters)[,paste0("cluster_",i)] <- array[,2])
    suppressWarnings(S4Vectors::mcols(clusters)[,paste0("DicerCall_",i)] <- array[,12])
    suppressWarnings(S4Vectors::mcols(clusters)[,paste0("Count_",i)] <- array[,4])
    suppressWarnings(S4Vectors::mcols(clusters)[,paste0("RPM_",i)] <- array[,5])
    # un-note line to inclde the major RNA strand.
    #suppressWarnings(S4Vectors::mcols(clusters)[,paste0("MajorRNA_",i)] <- array[,9])
  }
  data <- Repitools::annoGR2DF(clusters)
  data <- data %>% dplyr::select(-score, -phase, -source, -type)
  clustercolID <- data %>% dplyr::select(tidyselect::starts_with("cluster_"))
  clustercolID <- clustercolID[1]
  colnames(clustercolID)[1] <- 'clusterID'
  res <- data.frame(clusterID = clustercolID, data)

  if(!is.null(features)){
    overlap_data <- .import_annotations(features)
    res2 <- merge(data,overlap_data, by=c("chr","start", "end"),all.x=TRUE)
    clustercolID <- res2 %>% dplyr::select(tidyselect::starts_with("cluster_"))
    clustercolID <- clustercolID[1]
    colnames(clustercolID)[1] <- 'clusterID'
    res2 <- data.frame(clusterID = clustercolID, res2)
    return(res2)

  } else

    return(res)
}

