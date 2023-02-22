#' Quality control step for .fastq files
#'
#' @description Invokes an OS command to run unzip any compressed .gz files
#' and run a quality control analysis using \code{fastqc} on the files within the
#' directory.
#' @param  files The local directory to sRNA-seq sample. These can be raw or
#' trimmed.
#' @return A "qc" folder within the directory containing \code{fastqc} quality
#' check files for each sample in \code{HTML} format.
#' @examples
#' \dontrun{
#' raw_samples <- "/Users/user1/projectname/workplace/raw/"
#' qc <- checkQuality(raw_samples)
#' }
#' @export
#' @importFrom dplyr "%>%"
checkQuality <- function(files){
  dir.create(file.path(files,"qc"), showWarnings = FALSE)
  system(paste0("for i in ", files,"*.gz ; do gunzip $i; done"))
  system(paste0("for j in ",files,"*.fq; do fastqc $j -o ", files/qc/ "; done"))
}


