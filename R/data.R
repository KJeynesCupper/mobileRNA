#' The `sRNA_data` data set
#'
#' This data sets represents the results of the pre-processing. The data has
#' been organised using the `sample_table` function.
#'
#' Tomato (Solanium lycopersicum) and eggplant (Solanium melongena) were used
#' as grafting partners in a hetero-grafting experiment. The scion was
#' represented by tomato and the rootstock was represented by the eggplant. The
#' experiment used a tomato self-graft as a control.
#'
#' The dataset holds three replicates for the hetero-graft treatment (TomEgg_1,
#' TomEgg_2, TomEgg_3) and three replicates for the self-graft control
#' (TomTom_1, TomTom_2, TomTom_3).
#'
#' sRNA_data
#' @docType data
#' @keywords dataset
#' @name sRNA_data
#' @usage data(sRNA_data)
#' @description Simulated sRNAseq dataset
#' @details Simulates data is taken from eggplant and tomato sRNAseq samples and
#' created to simulate to movement of sRNA molecules from an Tomato rootstock to
#' an Eggplant Scion. Three Eggplant replicates were spiked with the same 150
#' tomato sRNA clusters, and named "heterograft_" 1 to 3. The analysis compares
#' these heterografts to three Eggplant self-grafts which are the original
#'  un-spiked Eggplant replicates.
#'  @examples
#'  data("sRNA_data")
NULL


#' sRNA_data_consensus
#' @docType data
#' @keywords dataset
#' @name sRNA_data_consensus
#' @usage data(sRNA_data_consensus)
#' @description Simulated sRNAseq dataset with RNA consensus
#' @details Simulates data is taken from eggplant and tomato sRNAseq samples and
#' created to simulate to movement of sRNA molecules from an Tomato rootstock to
#' an Eggplant Scion. Three Eggplant replicates were spiked with the same 150
#' tomato sRNA clusters, and named "heterograft_" 1 to 3. The analysis compares
#' these heterografts to three Eggplant self-grafts which are the original
#' un-spiked Eggplant replicates.
#' @examples
#'  data("sRNA_data_consensus")
NULL


#' sRNA_data_mobile
#' @docType data
#' @keywords dataset
#' @name sRNA_data_mobile
#' @usage data(sRNA_data_mobile)
#' @description Simulated sRNAseq dataset - potentially mobile RNAs
#' @details Simulates data is taken from eggplant and tomato sRNAseq samples and
#' created to simulate to movement of sRNA molecules from an Tomato rootstock to
#' an Eggplant Scion. Three Eggplant replicates were spiked with the same 150
#' tomato sRNA clusters, and named "heterograft_" 1 to 3. The analysis compares
#' these heterografts to three Eggplant self-grafts which are the original
#' un-spiked Eggplant replicates.
#' @examples
#'  data("sRNA_data_mobile")
NULL


#' chr12_Eggplant_V4.1.fa.gz
#' @docType data
#' @keywords dataset
#' @name chr12_Eggplant_V4.1.fa.gz
#' @description Simulated sRNAseq dataset - potentially mobile RNAs
#' @details Simulates data is taken from eggplant and tomato sRNAseq samples and
#' created to simulate to movement of sRNA molecules from an Tomato rootstock to
#' an Eggplant Scion. Three Eggplant replicates were spiked with the same 150
#' tomato sRNA clusters, and named "heterograft_" 1 to 3. The analysis compares
#' these heterografts to three Eggplant self-grafts which are the original
#' un-spiked Eggplant replicates.
#' @examples
#' system.file("extdata","chr12_Eggplant_V4.1.fa.gz", package="mobileRNA")
NULL

#' chr12_Eggplant_V4.1.fa.gz
#' @docType data
#' @keywords dataset
#' @format FASTA gzip
#' @name chr12_Eggplant_V4.1.fa.gz
#' @source Barchi, L., Rabanus‐Wallace, M. T., Prohens, J., Toppino, L.,
#' Padmarasu, S., Portis, E., ... & Giuliano, G. (2021). Improved genome
#' assembly and pan‐genome provide key insights into eggplant domestication and
#' breeding. The Plant Journal, 107(2), 579-596.
#' @references https://solgenomics.net/organism/Solanum_melongena/genome
#' @description Chromosome 12 of genome assembly Solanium melongena Version 4.1
#' @examples
#' system.file("extdata",
#' "chr12_Eggplant_V4.1.fa.gz", package="mobileRNA")
NULL

#' chr2_S_lycopersicum_chromosomes.4.00.fa.gz
#' @docType data
#' @keywords dataset
#' @format FASTA gzip
#' @name chr2_S_lycopersicum_chromosomes.4.00.fa.gz
#' @source Hosmani, P. S., Flores-Gonzalez, M., van de Geest, H., Maumus, F., Bakker, L.
#' V., Schijlen, E., ... & Saha, S. (2019). An improved de novo assembly and
#' annotation of the tomato reference genome using single-molecule sequencing,
#' Hi-C proximity ligation and optical maps. biorxiv, 767764.
#' @references https://solgenomics.net/organism/Solanum_lycopersicum/genome
#' @description Chromosome 2 of genome assembly Solanium lycopersicum (Heinz)
#' Version SL4.0
#' @examples
#' system.file("extdata",
#' "chr2_S_lycopersicum_chromosomes.4.00.fa.gz", package="mobileRNA")
NULL


#' chr2_ITAG4.0_gene_models.gff.gz
#' @docType data
#' @keywords dataset
#' @format GFF gzip
#' @name chr2_ITAG4.0_gene_models.gff.gz
#' @references https://solgenomics.net/organism/Solanum_lycopersicum/genome
#' @description Chromosome 2 of genome annotation Solanium lycopersicum (Heinz)
#' Version SL4.0
#' @examples
#' system.file("extdata",
#' "chr2_ITAG4.0_gene_models.gff.gz", package="mobileRNA")
NULL

#' chr12_Eggplant_V4.1_function_IPR_final.gff.gz
#' @docType data
#' @keywords dataset
#' @format GFF gzip
#' @source Barchi, L., Rabanus‐Wallace, M. T., Prohens, J., Toppino, L., Padmarasu, S.,
#' Portis, E., ... & Giuliano, G. (2021). Improved genome assembly and pan‐genome
#' provide key insights into eggplant domestication and breeding.
#' The Plant Journal, 107(2), 579-596.
#' @references https://solgenomics.net/organism/Solanum_melongena/genome
#' @name chr12_Eggplant_V4.1_function_IPR_final.gff.gz
#' @description Chromosome 12 of genome annotation Solanium melongena
#' Version 4.1
#' @examples
#' system.file("extdata",
#' "chr12_Eggplant_V4.1_function_IPR_final.gff.gz", package="mobileRNA")
NULL

