#' Create a working directory file structure
#'
#' @description Creates the basic directory to store your raw data,
#' genome references and annotation files
#'
#'
#' @param wd Path to working directory.
#'
#'
#' @return A working directory with folders the folders `"raw"` to store raw
#' sequencing data (`.fastq/.fq`), `"reference"` to store reference
#' genomes (`.fa/.fasta`) & `"annotation"` to store reference
#' annotation files (`.gff/.gff3`)
#'
#' @details This argument could include using the [base::getwd()] function to
#' retrieve path to the working directory.
#'
#' @examples
#' \dontrun{
#' createWorkplace(getwd())
#' }
#' @export
createWorkplace <- function(wd){
  message("The working directory path is ", wd)
  message("The following folders have been created here: reference,
          annotation, raw")
  dir.create(file.path(wd, "/reference"), showWarnings = FALSE)
  dir.create(file.path(wd, "/reference/merge"), showWarnings = FALSE)
  dir.create(file.path(wd, "/annotation"), showWarnings = FALSE)
  dir.create(file.path(wd, "/annotation/merge"), showWarnings = FALSE)
  dir.create(file.path(wd, "/raw"), showWarnings = FALSE)
}

