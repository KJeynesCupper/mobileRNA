#' Sample Data Table
#'
#' @description Import specific data from the \code{"Results"} file of each
#' sample and collate into a single dataframe. Plus, calculates the \code{FPKM}
#'  for each sample.
#' @param results Local directory location of sample folders. NOTE: Following
#' the pre-processing steps, these can be found in second alignment folders.
#' It is recommended to use results from the unique alignment found in
#' \code{"./second_alignment_unique/"} over multimapped reads.
#' It is recommended to run a parallel analysis one using the uniquely
#' mapped reads and the other using multimapped reads.
#' @param samples Vector of characters naming the sample names correlating
#' to outputted folders either \code{"second_alignment_unique"} or
#' \code{"second_alignment_multi"}.
#' @param clusters Import \code{.gff3} from \code{"reference"} folder.
#' @param totalNumReads Table containing containing the total number of reads
#' for each samples. Each value labeled with the sample name.
#'
#' @param features Path to a .bed file containing the genomic information
#' (annotations/features) which overlaps with the locations of the
#' dicer-derived clusters. The file is created as a part of the pre-processing
#' steps, see fun[find_overlap] function for more information.
#' @return A large dataframe with data for each sample including Dicercall,
#' Count and Cluster ID.
#' @examples
#' \dontrun{
#' genomic_features <- "./data/reference/overlap.bed"
#' directory <- "./data/alignment_unique_two/"
#' sample_names <- c("TomEgg_1","TomEgg_2", "TomEgg_3",
#'               "TomTom_1", "TomTom_2", "TomTom_3" )
#'
#' clusters <- rtracklayer::import.gff("./data/reference/ClustersInfo.gff3")
#'
#' total_reads <- c(24441759,21378845, 21482356, 3951725, 3343954, 2586910)
#' names(total_reads) <- samples
#' siRNA_data <- RNAimport(results = directory,
#'                            samples = sample_names,
#'                            clusters = clusters,
#'                            totalNumReads = total_reads
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

RNAimport <- function(results, samples, clusters,
                         totalNumReads, features){
  if (base::missing(results)|| !base::inherits(results, c("character"))) {
    stop("results must be an object of character, representing a path to a
         directory")
  }
  if (base::missing(samples)|| !base::inherits(samples, c("character"))) {
    stop("samples must be an object of character, representing a vector of
         sample names")
  }
  #if (base::missing(clusters)|| !base::inherits(clusters, c("character"))) {
    #stop("clusters must be an object of character, representing a path to a
      #   directory") }
  if (base::missing(totalNumReads)||
      !base::inherits(totalNumReads,c("numeric"))) {
    stop("totalNumReads must be an object of named numeric")
  }
  if (missing(features)|| !base::inherits(features, c("character"))) {
    stop("features must be an object of character, representing a path to a
         directory")
  }
  for (i in samples) {
    array <- utils::read.table(paste0(results,i,"/Results.txt", sep=""),
                               sep="\t", header =F, row.names=NULL)
    S4Vectors::mcols(clusters)[,paste0("DicerCall_",i)] <- array[,12]
    S4Vectors::mcols(clusters)[,paste0("Count_",i)] <- array[,4]
    S4Vectors::mcols(clusters)[,paste0("cluster_",i)] <- array[,2]
  }
 for (j in samples) {
   S4Vectors::mcols(clusters)[,paste0("FPKM_",j)] <-
     S4Vectors::mcols(clusters)[,paste0("Count_",j)] /
     (BiocGenerics::width(clusters)/1000 * totalNumReads[j]/1000000 )
 }
  data <- Repitools::annoGR2DF(clusters)
  overlap_data <- .import_annotations(features)
  data$features <- overlap_data$features[base::match(data$start, data$start)]
  data <- data %>% dplyr::select(-score, -phase, -source, -type)
  return(data)
}



