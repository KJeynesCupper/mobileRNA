#' Make the working directory structure
#'
#' @description Creates the basic directory structure to store your raw data,
#' genome references and annotation files
#' @param  dirwd Working Directory path. This argument could include using
#' the [getwd()] function to retrieve path.
#' @return A working directory with folders the folders `"raw"` to store raw
#' sequencing data (`.fastq/.fq`), `"reference"` to store reference
#' genomes (`.fa/.fasta`) & `"annotation"` to store reference
#' annotation files (`.gff/.gff3`)
#' @examples
#' \dontrun{
#' create_workplace(getwd())
#' }
#' @export
create_workplace <- function(dirwd){
  message("The working directory path is ", getwd())
  message("The following folders have been created here:reference,
          annotation, raw", getwd())
  dir.create(file.path(dirwd, "/reference"), showWarnings = FALSE)
  dir.create(file.path(dirwd, "/reference/merge"), showWarnings = FALSE)
  dir.create(file.path(dirwd, "/annotation"), showWarnings = FALSE)
  dir.create(file.path(dirwd, "/annotation/merge"), showWarnings = FALSE)
  dir.create(file.path(dirwd, "/raw"), showWarnings = FALSE)
}

