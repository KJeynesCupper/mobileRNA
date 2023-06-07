#' Merge two FASTA reference genomes into one
#'
#' @description Merges two reference genomes (.fa/.fasta). into one single
#' reference with modified chromosome names to ensure distinguishability.
#'
#'
#'@param genomeA the path  \code{base::connection()} to a genome reference
#'file in FASTA format (.fa/.fasta). File may be supplied as compressed with a
#'.gzip extension.
#'
#'@param genomeB the path  \code{base::connection()} to a genome reference
#'file in FASTA format (.fa/.fasta). File may be supplied as compressed with a
#'.gzip extension.
#'
#'@param out_dir either a character string or a \code{base::connections()} open
#'for writing. Place path to output directory in "", including file output name.
#'
#'
#'@param abbreviationGenomeA a string placed in "", to replace chromosome names
#'within \code{genomeA}.Default set as "A".
#'
#'@param abbreviationGenomeB a string placed in "", to replace chromosome names
#'within \code{genomeB}.Default set as "B".
#'
#'@param replace_chr_names a logical value indicating whether the chromosome
#'names of the supplied annotation files are to be altered or not. As default,
#'chromosome names are altered.
#'
#'
#'@return A FASTA format file containing two genome assemblies with edited
#'chromosome names (prefixes, and removal of periods)
#'@details
#'
#' The functions primary goal is to merge two FASTA files, however, when
#' merging genomic files it is critical that the two genomes are distinguishable
#' by the chromosome names. As a default setting, the function extracts the
#' chromosome names for a given FASTA file and alters the name while retaining
#' the identifying number.
#'
#' The function requires the input of two reference genomes, where one
#' represents `Genome-A` and the other represents `Genome-B`. As default, the
#' function will rename the chromosome names in `Genome-A` to "A". For example,
#' A0, A1, A2 etc. To set a custom chromosome name for `Genome-A` alter the
#' argument \code{abbreviationGenomeA}. While, for  `Genome-B` as default the
#' chromosome names will be named "B", for example, B0, B1, B2 etc. To set a
#' custom chromosome name for `Genome-A` alter the argument
#' \code{abbreviationGenomeB}. The function can do so by  draw the chromosome
#' number within the given GFF file, remove all prior character or numerical
#' values, and replace it with the supplied string.
#'
#'
#'IMPORTANT:  The genome reference and annotation of a species
#'must have chromosomes with matching names. It is critical that if you used
#'the [mobileRNA::RNAmergeAnnotations()] function to to create a merged genome
#'annotation,that you treat the input references in the same way.
#'
#' @examples \dontrun{
#' merged_ref <- RNAmergeGenomes(genomeA = "./workplace/reference/ref1.fa",
#'             genomeB = "./workplace/reference/ref2.fa",
#'             out_dir = "./workplace/reference/merge/merged_ref.fa")
#'
#' ## or, to set specific changes to chromosome names. annotationA represents
#' ## the Solanum lycopersicum and the chromosomes will be abbreviated to `SL`,
#' ## and annotationB represents Solanum melongena and the chromosomes will be
#' ## abbreviated to `SM`.
#'
#' merged_ref2 <- RNAmergeGenomes(genomeA = "./workplace/reference/ref1.fa",
#'             genomeB = "./workplace/reference/ref2.fa",
#'             out_dir = "./workplace/reference/merge/merged_ref.fa",
#'             abbreviationGenomeA = "SL",
#'             abbreviationGenomeB = "SM")
#'
#'}
#'
#'
#' @importFrom Biostrings "readDNAStringSet"
#' @importFrom tools "file_ext"
#' @importFrom Biostrings "writeXStringSet"
#'
#' @export
# Load required packages
RNAmergeGenomes <- function(genomeA, genomeB,
                            out_dir,
                            abbreviationGenomeA= "A",
                            abbreviationGenomeB= "B",
                           replace_chr_names = TRUE) {

  if (base::missing(genomeA) || !base::inherits(genomeA, c("character"))) {
    stop(paste("Please specify genomeA, a connection to a FASTA file in local"))
  }

  if (base::missing(genomeB) || !base::inherits(genomeB, c("character"))) {
    stop(paste("Please specify annotationA, a connection to a FASTA file in local"))
  }

  if (base::missing(out_dir) || !base::inherits(out_dir, c("character")) || tools::file_ext(out_dir == c("fa", "fasta"))) {
    stop(paste("Please specify out_dir, a connection to a local directory to
               write and save merged annotation. Ensure file name with extension
               (.fa or .fasta) is supplied."))
  }


  message("Loading reference genomes")

  # Read reference genomes
  ref1 <- Biostrings::readDNAStringSet(file.path(genomeA), format = "fasta")
  ref2 <- Biostrings::readDNAStringSet(file.path(genomeB), format = "fasta")

  # Replace chromosome names in reference genomes


  if (replace_chr_names) {
    message("Replacing chromosome names")
    ref1_names <- names(ref1)
    ref1_newnames <- paste0(abbreviationGenomeA,"_", ref1_names )
    ref1_newnames <- sub("\\.", "", ref1_newnames)
    names(ref1) <- ref1_newnames

    ref2_names <- names(ref1)
    ref2_newnames <- paste0(abbreviationGenomeB,"_", ref1_names )
    ref2_newnames <- sub("\\.", "", ref2_newnames)
    names(ref2) <- ref2_newnames
  }

  ref1_save <- paste0(tools::file_path_sans_ext(genomeA), "_altered.fa")
  ref2_save <- paste0(tools::file_path_sans_ext(genomeB), "_altered.fa")


  message("Writting altered reference files")

  Biostrings::writeXStringSet(ref1,ref1_save , format = "fasta", append = FALSE)
  Biostrings::writeXStringSet(ref2, ref2_save, format = "fasta", append = TRUE)


  message("Writting merged reference file to:", out_dir)
  system(paste0("cat ", ref1_save, " ", ref2_save, " >", out_dir))

  }



