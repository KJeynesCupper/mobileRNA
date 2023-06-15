#'Extract sRNA cluster sequences
#'
#'
#'
#'@description Extract sRNA cluster sequences; identifies whether the most
#'abundant sRNA for a sRNA cluster id consistent across the replicates, and if
#'so extracts the sRNA nucleotide sequence.
#'
#' The function also calculates the RNA & DNA complementary sequences, as well
#' as stating the length of the sequence.
#'
#' The results can be stored in the the working dataframe supplied to the
#' function or stored as a new dataframe containing only the results from the
#' function. This includes:
#' - Match:  logical; whether the sRNA sequence is consistent across replicates
#' - Sequence:  character; sequence of the most abundant sRNA within a cluster
#' - Complementary_RNA: character; complementary RNA nucleotide sequence
#' - Complementary_DNA: character; complementary DNA nucleotide sequence
#' - Width: numeric; length of nucleotide sequence
#'
#'
#'@param data dataframe; containing columns starting with `Cluster` & `MajorRNA`
#'
#'@param original logical; output results to working date as additional
#'columns (`original=TRUE`), or output as new dataframe (`original=FALSE`).
#'by default, FALSE
#'
#'
#'@return The results can be stored in the the working dataframe supplied to the
#' function or stored as a new dataframe containing only the results from the
#' function. The results includes:
#' - Match:  logical; whether the sRNA sequence is consistent across replicates
#' - Sequence:  character; sequence of the most abundant sRNA within a cluster
#' - Complementary_RNA: character; complementary RNA nucleotide sequence
#' - Complementary_DNA: character; complementary DNA nucleotide sequence
#' - Width: numeric; length of nucleotide sequence
#'
#'
#'
#'@examples
#'
#'data('sRNA_data_mobile')
#'
#'mobile_sequences <- RNAsequences(sRNA_data_mobile)
#'
#'@importFrom dplyr "select"
#'@importFrom dplyr "starts_with"
#'@importFrom dplyr "where"
#'@importFrom dplyr "mutate"
#'@export


RNAsequences <- function(data, original = FALSE){

  # select only columns with RNA seqs, remove columns with only NA values
  df <- data %>%
    dplyr::select(dplyr::starts_with("MajorRNA")) %>%
    dplyr::select(-dplyr::where(~all(. == "N")))
  # add cluster names


  # add new columns, if match then true or false, then add sequence consensus.
  df$Match <- apply(df, 1, function(row) all(row == row[1]))
  df$Sequence <- ifelse(df$Match, df[,2], NA)
  df$width <- nchar(df$Sequence)

  # calculate complementary sequences
  df$Complementary_RNA <- sapply(df$Sequence, find_complementary_sequenceRNA)
  df$Complementary_DNA <- sapply(df$Sequence, find_complementary_sequenceDNA)
  # add as col
  Cluster <-  data$Cluster
  df <- cbind(Cluster, df )

  # add columns to og data
  if (original){
    data_output <- data %>%
      dplyr::mutate(Match = df$Match,
             Sequence = df$Sequence,
             Width = df$width,
             Complementary_RNA = df$Complementary_RNA,
             Complementary_DNA = df$Complementary_DNA)
  }else # made new df
    if(original == FALSE){
      data_output <- df %>%
        dplyr::select(!dplyr::starts_with("MajorRNA"))
    }
  return(data_output)
}





