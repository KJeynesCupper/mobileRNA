#' Quality control of FASTQ files
#'
#' @description Invokes an OS command to \code{fastqc} to run quality control
#' analysis on the files within the choosen directory.
#'
#' @details The analysis can have compressed or uncompressed FASTQ files as
#' input. If files are compressed in the \code{.gzip} format, these will be
#' unzipped, and a copy will be retained in the directory.
#'
#' This can be utilised before and after the trimming of samples. Checking the
#' quality of sequenced samples before analysis will identify whether adapter
#' sequences or low quality reads need to be removed to improve the sample
#' quality in preparation for the analysis. Similarly, checking sample quality
#' after trimming will ensure the trimming step was performed as intended.
#'
#'
#' @param files The directory to FASTQ samples.
#'
#' @return Outputs two files per sample, one in html format and the other in
#' json format. These are stored in a new folder within the directory containing
#' the FASTQ files; this folder is names "qc".
#'
#' @examples
#' \dontrun{
#' raw_samples <- "/Users/user1/projectname/workplace/raw/"
#' qc <- checkQuality(raw_samples)
#' }
#' @export
checkQuality <- function(files){
  dir.create(file.path(files,"qc"), showWarnings = FALSE)
  system(paste0("for i in ", files,"*.gz ; do gunzip $i; done"))
  system(paste0("for j in ",files,"*.fq; do fastqc $j -o ", files/qc/ "; done"))
}




