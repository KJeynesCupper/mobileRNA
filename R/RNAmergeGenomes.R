#' Merge two FASTA genome assemblies
#'
#' @description Merges two reference genomes (.fa/.fasta). into one single
#' reference with modified chromosome names to ensure distinguishability.
#'
#'
#'@param genomeA a large DNAStringSet; a genome reference assembly file in
#'FASTA format (.fa/.fasta).
#'
#'@param genomeB a large DNAStringSet; a genome reference assembly file in
#'FASTA format (.fa/.fasta).
#'
#'@param out_dir either a character string or a \code{base::connections()} open
#'for writing. Place path to output directory in "", including file output name.
#'Output name must have a ".fa" extension.
#'
#'@param abbreviationGenomeA a string placed in "", to replace chromosome names
#'within \code{genomeA}.Default set as "A".
#'
#'@param abbreviationGenomeB a string placed in "", to replace chromosome names
#'within \code{genomeB}.Default set as "B".
#'
#'
#'@return Returns a single FASTA format file containing both  genome assemblies
#'with edited chromosome names (prefixes, and removal of periods) to the give
#'directory, as well as the individual edited genomes assemblies named
#'"genomeA_altered.fa" and "genomeB_altered.fa", respectively.
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
#' @examples
#'
#'# import FASTA files into R
#'
#'  # load into R
#'  ref1 <- Biostrings::readDNAStringSet(system.file("extdata",
#'                            "chr12_Eggplant_V4.1.fa.gz", package="mobileRNA"),
#'                             format = "fasta")
#'
#'  ref2 <- Biostrings::readDNAStringSet(system.file("extdata",
#'     "chr2_S_lycopersicum_chromosomes.4.00.fa.gz", package="mobileRNA"),
#'      format = "fasta")
#'
#' # run function to merge
#' merged_ref <- RNAmergeGenomes(genomeA = ref1, genomeB = ref2,
#'                                out_dir = "../references/merged_ref.fa")
#'
#'
#'
#'
#' ## or, to set specific changes to chromosome names:
#'  # genomeA represents Solanum melongena and the chromosomes will be
#'  # abbreviated to `SM.
#' # genomeB represents Solanum lycopersicum and the chromosomes will be
#' # abbreviated to `SL`.
#'
#' merged_ref2 <- RNAmergeGenomes(genomeA = ref1,
#'             genomeB = ref2 ,
#'             out_dir = "../references/merged_ref.fa",
#'             abbreviationGenomeA = "SM",
#'             abbreviationGenomeB = "SL")
#'
#' @importFrom Biostrings "readDNAStringSet"
#' @importFrom Biostrings "writeXStringSet"
#'
#' @export
# Load required packages
RNAmergeGenomes <- function(genomeA, genomeB,
                            out_dir,
                            abbreviationGenomeA= "A",
                            abbreviationGenomeB= "B") {

  if (base::missing(genomeA)) {
    stop(paste("Please specify genomeA, a connection to a FASTA file in local"))
  }

  if (base::missing(genomeB)) {
    stop(paste("Please specify annotationA, a connection to a FASTA file in local"))
  }

  if (base::missing(out_dir) || !grepl("\\.fa$", out_dir)) {
    stop(paste("Please specify out_dir, a connection to a local directory to
               write and save merged annotation. Ensure file name with extension
               (.fa or .fasta) is supplied."))
  }

  # Replace chromosome names in reference genomes
    message("Replacing chromosome names")
    ref1_names <- names(ref1)
    ref1_newnames <- paste0(abbreviationGenomeA,"_", ref1_names )
    ref1_newnames <- sub("\\.", "", ref1_newnames)
    names(ref1) <- ref1_newnames

    ref2_names <- names(ref2)
    ref2_newnames <- paste0(abbreviationGenomeB,"_", ref2_names )
    ref2_newnames <- sub("\\.", "", ref2_newnames)
    names(ref2) <- ref2_newnames

# location to save
  location <- dirname(out_dir)
  ref1_save <- paste0(location,"/", "genomeA_altered.fa")
  ref2_save <- paste0(location,"/", "genomeB_altered.fa")

  message("Writting altered reference files to output location: ", location)

  Biostrings::writeXStringSet(ref1,ref1_save , format = "fasta", append = FALSE)
  Biostrings::writeXStringSet(ref2, ref2_save, format = "fasta", append = TRUE)


  message("Writting merged reference file to: ", out_dir)
  system(paste0("cat ", ref1_save, " ", ref2_save, " >", out_dir))

  }






















