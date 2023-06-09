#' Distribution of sRNA cluster loci across genomic features
#'
#' @description Calcularepeats the percentage of sRNA dicer-derived clusters
#' which overlap with genomic features, including 2-kb promoter regions, repeat
#' regions, exons, introns, and 3'/'5 untranslated regions.
#'
#'
#' @details
#' The function calculates the number or percentage of sRNA dicer-derived
#' clusters which overlap with genomic features based on their genomic
#' coordinates. This is outputted as either a matrix or data frame,
#' respectively.
#'
#' This function can be utilised at different steps in your analysis, but is
#' particularly powerful when observing the genomic location of potential mobile
#' sRNA, and can be overlapped with a specific genome rather than a merged
#' annotation file.
#'
#' @param data data frame; containing rows of potential dicer-derived clusters
#' including columns which supply the genomic coordinates, where `chr` supplies
#' the chromosome number, `start` and `end` which supply the start and
#' end coordinates.
#'
#' @param annotation A path, URL, connection or GFFFile object. A genome
#' reference annotation file (.gff/.gff1/.gff2/.gff3).
#'
#' @param repeats A path, URL, connection or GFFFile object. A genome
#' reference annotation file, which specifically contains information on repeat
#' sequences in the genome (.gff/.gff1/.gff2/.gff3). By default, this is not
#' required, however if there is a specific repeats annotation file for the
#' genome it is suggested to supply it.
#'
#'
#' @param promoterRegions numeric; define promoter region upstream of genes.
#' Default is `promoterRegions=2000` , ie promoters set at 2kb upstream of genes
#' @param percentage returns results as a percentage of the total, when
#' \code{percentage = TRUE} (default). While \code{percentage = FALSE}, results
#' are returned as a count value representing the number of sRNA clusters
#' that overlap with a given genomic feature.
#'
#'
#'
#' @return Returns a table containing the number or percentage of overlaps in
#' the supplied dataset with specific regions in the genome annotation such
#' as genes, repeats, introns, exons.
#'
#' @examples
#' \dontrun{
#'
#' dis_features <- RNAfeatures(data = sRNA_data_consensus,
#'                        annotation = "./annotation/eggplant_genome.gff3",
#'                        repeats = "./annotation/eggplant_genome_repeats.gff3")
#' }
#'@importFrom rtracklayer "import.gff3"
#'@importFrom GenomicRanges "setdiff"
#'@importFrom stats "start"
#'@importFrom stats "end"
#'@importFrom BiocGenerics "strand"
#'@importFrom dplyr "select"
#'@importFrom dplyr "mutate"
#'@importFrom dplyr "filter"
#'@importFrom GenomicRanges "makeGRangesFromDataFrame"
#'@importFrom scales "label_percent"
#'@importFrom dplyr "%>%"
#'@importFrom IRanges "overlapsAny"
#'@importFrom BiocGenerics "width"
#'@importFrom Repitools "annoGR2DF"
#'@export
RNAfeatures <- function(data, annotation,
                        repeats = NULL,
                        promoterRegions = 2000,
                        percentage = TRUE){
  annotation_info <-rtracklayer::import.gff3(annotation)

  anno_repeats <- repeats

  if(is.null(anno_repeats)) {
    # features
    repeats <-subset(annotation_info, type=="transposable_element",
                type=="transposable_element_gene")
    } else {
      anno_repeats < rtracklayer::import(anno_repeats)
      repeats <-subset(anno_repeats, type=="transposable_element",
                  type=="transposable_element_gene")
    }
  genes<-subset(annotation_info[!IRanges::overlapsAny(annotation_info,
                                                            repeats,
                                                          ignore.strand=TRUE)],
                                                             type=="gene")
  five_UTR<-subset(annotation_info[IRanges::overlapsAny(annotation_info,
                                                              genes,
                                                          ignore.strand=TRUE)],
                                                         type=="five_prime_UTR")
  five_UTR<-GenomicRanges::setdiff(five_UTR, c(repeats), ignore.strand=TRUE)
  three_UTR<-subset(annotation_info[IRanges::overlapsAny(annotation_info,
                                                               genes,
                                                          ignore.strand=TRUE)],
                                                       type=="three_prime_UTR")
  three_UTR<-GenomicRanges::setdiff(three_UTR, c(five_UTR, repeats),
                                    ignore.strand=TRUE)
  exons<-subset(annotation_info[IRanges::overlapsAny(annotation_info,
                                                           genes,
                                                           ignore.strand=TRUE)],
                type=="exon")
  exons<-GenomicRanges::setdiff(exons, c(five_UTR, three_UTR, repeats),
                                ignore.strand=TRUE)
  introns <- GenomicRanges::setdiff(genes, c(exons, five_UTR, three_UTR),
                                    ignore.strand=TRUE)

  # define promoter regions
  gene_promoters <-Repitools::annoGR2DF(genes)
  pos_strand_promoter <- gene_promoters %>%
    dplyr::filter(strand == "+") %>% dplyr::mutate(end=start) %>%
    dplyr::filter(strand == "+") %>% dplyr::mutate(start=start-promoterRegions)

  neg_strand_promoter <- gene_promoters %>%
    dplyr::filter(strand == "-") %>% dplyr::mutate(end=start) %>%
    dplyr::filter(strand == "-") %>% dplyr::mutate(start=start-promoterRegions)

  promoters <- rbind(pos_strand_promoter, neg_strand_promoter)
  promoters <- GenomicRanges::makeGRangesFromDataFrame(promoters)


  promoters <-   BiocGenerics::setdiff(promoters, c(repeats, exons,five_UTR,
                                                    three_UTR,introns),
                                       ignore.strand=TRUE)
  #get others
  others <-   BiocGenerics::setdiff(subset(annotation, type=="chromosome"),
                                    c(repeats, exons,five_UTR,three_UTR,introns,
                                      promoters), ignore.strand=TRUE)
  # data frame
  sRNA_features_df <- matrix(0, ncol=7, nrow = 2)
  colnames(sRNA_features_df) <- c("promoters","exons", "introns", "5'UTR",
                                  "3'UTR", "repeats", "others")
  rownames(sRNA_features_df) <- c("Genome", "Dataset")

  # genome
  sRNA_features_df[1,1] <- sum(BiocGenerics::width(promoters))
  sRNA_features_df[1,2] <- sum(BiocGenerics::width(exons))
  sRNA_features_df[1,3] <- sum(BiocGenerics::width(introns))
  sRNA_features_df[1,4] <- sum(BiocGenerics::width(five_UTR))
  sRNA_features_df[1,5] <- sum(BiocGenerics::width(three_UTR))
  sRNA_features_df[1,6] <- sum(BiocGenerics::width(repeats))
  sRNA_features_df[1,7] <- sum(BiocGenerics::width(others))

  # select sample
  sRNA_df <-  data %>% dplyr::select(chr, start, end)
  sRNA_df <-  GenomicRanges::makeGRangesFromDataFrame(sRNA_df)

  sRNA_features_df[2,1] <- sum(BiocGenerics::width(intersect(
    promoters,sRNA_df)))
  sRNA_features_df[2,2] <- sum(BiocGenerics::width(intersect(
    exons,sRNA_df)))
  sRNA_features_df[2,3] <- sum(BiocGenerics::width(intersect(
    introns,sRNA_df)))
  sRNA_features_df[2,4] <- sum(BiocGenerics::width(intersect(
    five_UTR,sRNA_df)))
  sRNA_features_df[2,5] <- sum(BiocGenerics::width(intersect(
    three_UTR,sRNA_df)))
  sRNA_features_df[2,6] <- sum(BiocGenerics::width(intersect(repeats,sRNA_df)))
  sRNA_features_df[2,7] <- sum(BiocGenerics::width(intersect(others,sRNA_df)))

  if(percentage == TRUE){
    # convert to percentage
    sRNA_features_df <- data.frame(t(sRNA_features_df)) %>%
      dplyr::mutate(Genome = scales::label_percent()(Genome / sum(Genome)))%>%
      dplyr::mutate(Dataset = scales::label_percent()(Dataset / sum(Dataset)))
    return(sRNA_features_df)
  } else
  return(sRNA_features_df)
}
