#' Identify potential mobile sRNA molecules
#'
#' @description A function to identify potential mobile small RNA molecules
#' traveling from one genotype to another in a hetero-grafted system.
#'
#' @details
#' The function undertakes two different roles. First, it selects the sRNA
#' cluster which are mapped to a particular genome. It does so by either keeping
#' or removing sRNA mapped to chromosomes. Hence, this function will only work
#' if the two genomes are distinguishable by their chromosome names.
#'
#' The second step, after selecting clusters which are mapped to the genome of
#' interest, removes sRNA clusters which were incorrectly mapped. These
#' are clusters which have counts or FPKM values in the control samples
#' (ie. same genome as the destination tissue in the hetero-graft condition).
#' These samples should not have counts if the sRNA originates from a different
#' genotype to the control.
#'
#' The function also allows for statistical analysis based on the results
#' collect from differential analysis of the total dataset using the function
#' [RNAlocate::RNAanalysis()]. This features enables the filtering of sRNA
#' clusters which meet a specific p-value or adjusted p-values.
#'
#' @param data Numeric data frame
#'
#' @param controls Character vector; containing names of control samples.
#'
#' @param id a string related to the chromosomes in a particular genome. A
#' distinuishing feature of the genome of interest or non-interest in the
#' chromosome name (`chr` column).
#'
#' @param task an option to keep or remove the chromosomes containing the
#' identifying string. To keep the chromosomes with the ID, set task=keep.
#' To remove, set `task="remove"`.
#'
#'
#' @param statistical If TRUE, will undertake statistical filtering based on the a
#' p-value or adjusted p-value threshold stated by `padj` and `p.value`.
#'
#' @param padj A user defined numeric value to represent the adjusted p-value
#' threshold to define statistic significance. Defaults set at 0.05.Only mobile
#' molecules with adjusted p-values equal or lower than specified are returned.
#'
#' @param p.value A user defined numeric value to represent the p-value
#' threshold to define statistic significance. There is no default value, set
#' this instead of using an adjusted p-value to filter molecules. Only mobile
#' molecules with p-values equal or lower than specified are returned.
#'
#'
#' @return A data-frame containing statistically significant mobile small
#' RNA molecules.
#'
#' @examples
#'
#'
#' ## DESeq2 example (24-nt & 21/22-nt sRNA)
#'data("sRNA_24_DESeq2")
#'data("sRNA_1224_DESeq2")
#'
#'
#' # vector of control names
#' controls <- c("TomTom_1", "TomTom_2", "TomTom_3")
#'
#' # Identify mobile 24-nt sRNA from eggplant genome:
#' # remove clusters associated to tomato
#' sRNA_24_mobile_DEseq2 <- RNAmobile(data = sRNA_24_DESeq2,
#'                                          controls = controls,
#'                                          id = "SL40",
#'                                          task = "remove")
#'
#' # Identify mobile 21/22-nt sRNA from eggplant genome:
#' # remove clusters associated to tomato
#' sRNA_2122_mobile_DEseq2  <- RNAmobile(data = sRNA_2122_DESeq2,
#'                                              controls = controls,
#'                                              id = "SL40",
#'                                              task = "remove")
#'
#'
#'
#'
#'
#'
#' ## Mobile molecules from the edgeR example
#'
#' data("sRNA_24_edgeR")
#' data("sRNA_2122_edgeR")
#' # vector of control names
#' controls <- c("TomTom_1", "TomTom_2", "TomTom_3")
#'
#'
#' # Identify mobile 24-nt sRNA from eggplant genome:
#' # remove clusters associated to tomato
#' sRNA_24_mobile_edgeR <- RNAmobile(data = sRNA_24_edgeR,
#'                                          controls = controls,
#'                                          id = "SL40",
#'                                          task = "remove")
#'
#' # Identify mobile 21/22-nt sRNA from eggplant genome:
#' # remove clusters associated to tomato
#' sRNA_2122_mobile_edgeR <- RNAmobile(data = sRNA_2122_edgeR,
#'                                            controls = controls,
#'                                            id = "SL40",
#'                                            task = "remove")

#'
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom dplyr "filter"
#' @importFrom dplyr "select"
#' @importFrom tidyselect "starts_with"

RNAmobile <- function(data,controls, id, task =c("keep", "remove"),
                             statistical = TRUE, padj = 0.05, p.value = NULL){
  if (base::missing(task) || !task %in% c("keep", "remove")) {
    stop(paste("Please specify task as to either", "(\"keep\", or \"remove\")",
               "chromosomes with the coresponding string"))
  }
  if (!base::inherits(data, c("matrix", "data.frame", "DataFrame"))) {
  stop("data must be an object of class matrix, data.frame, DataFrame")
    }
  if (base::missing(id) || id %in% "") {
    stop(paste("Please specify a single character string which is present in
               the all the chromosomes within the genome you wish to keep
               or remove"))
  }
  if (base::missing(id) || id %in% "") {
    stop(paste("Please specify a single character string which is present in
               the all the chromosomes within the genome you wish to keep
               or remove"))
  }
  if(task== 'remove'){
    x <- data %>% dplyr::filter(!base::grepl(id,chr))
  } else if (task == "keep"){
      x <- data %>% dplyr::filter(base::grepl(id,chr))
    }else
      stop(base::sQuote(x), " not implemented")
  y <- .remove_mapping_errors(data = x, controls = controls)
  mean_FPKM <- base::rowMeans(y %>%
                                dplyr::select(tidyselect::starts_with("FPKM")))
  mean_count <- base::rowMeans(y %>%
                                 dplyr::select(tidyselect::starts_with("Count")))
  CPM <- base::rowMeans(y %>% dplyr::select(tidyselect::starts_with("Count")))
  y$FPKM_mean <- mean_FPKM
  y$Count_mean <- mean_count
  if (statistical) {
    if (is.null(p.value)) {
      res <- y %>% filter(padjusted <= padj)
    } else
      res <- y %>% filter(pvalue <= p.value)
  } else
    if (statistical == FALSE){
      res <- y
    }
  return(res)
}

