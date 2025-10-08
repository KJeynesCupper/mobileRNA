#' Identify putative RNA molecules produced by the non-tissue sample genome
#'
#' @description A function to identify the putative sRNA or mRNA  molecules
#' produced by the non-tissue sample genome. Includes putative RNA mobilome or
#' RNAs not expected to be found within the tissue of origin.
#'
#'@references
#' baymobil \url{https://doi.org/10.21203/rs.3.rs-2520491/v1}
#'
#'
#' @details
#' The function identifies candidate sRNAs or mRNAs produced by a specific
#' genome/genotype. It does so by either keeping or removing those mapped to a
#' given genome. To do so, it requires a common pre-fix across chromosomes of
#' the given genome. See[mobileRNA::RNAmergeGenomes()] for more information.
#' In addition, it removes RNAs which were likely to be falsely mapped. These
#' are those which were mapped to the non-tissue genotype in the control
#' samples.
#'
#' **For sRNAseq:**
#' A greater confidence in the sRNA candidates can be achieved by setting
#' a threshold that considers the number of replicates which contributed to
#' defining the consensus dicercall (ie. consensus sRNA classification). This
#' parameter filters based on the `DicerCounts` column introduced by the
#' [mobileRNA::RNAdicercall()] function.
#'
#' **For mRNAseq:**
#' A greater confidence in the mRNA candidates can be achieved by setting
#' a threshold that considers the number of replicates which contained reads
#' for the mRNA molecule. This parameter filters based on the `SampleCounts`
#' column introduced by the [mobileRNA::RNAimport()] function.
#'
#' Alternatively, utilize the `baymobil` algorithm, which identifies the mobile
#' mRNAs using exact Bayesian inference based on SNP variant information. This
#' option requires the user to input substantially more information. Please
#' be aware this option requires significant processing time (so, maybe get a
#' brew while this runs). Note that this option produces additional files, see
#' `Value` section for more details.
#'
#' To utilise `baymobil` please ensure that you have installed baymobil and
#' added it to your path. As well as a conda environment containing `bcftools`
#' and `snpEff`. This is essential. Pleasue supply the conda environment
#' that includes these packages to the `condaenv` parameter.
#'
#'
#' **Statistical Analysis**
#' The function also allows for filtering using statistical inference generated
#' from the differential analysis of the total data set using the function
#' [mobileRNA::RNAdifferentialAnalysis()]. When `statistical=TRUE`, the feature
#' is enabled and selects molecules that meet the adjusted p-value
#' cutoff defined by `alpha`.
#'
#'
#'
#' @param input character; must be either "sRNA" or "mRNA" to represent the type
#' of data.
#'
#' @param data data.frame; generated through the \pkg{mobileRNA} method.
#'
#' @param controls character vector; containing names of control samples.
#'
#' @param genome.ID character; string or chromosome identifier related to the
#' chromosomes in a given genome. A distinguishing feature of the genome of
#' interest or non-interest in the chromosome name (`chr` column).
#'
#' @param task character; string to set the method to keep or remove the
#' chromosomes containing the identifying string. To keep the chromosomes with
#' the ID, set task=keep. To remove, set `task="remove"`. As default, task is
#' set to `keep`.
#'
#' @param statistical If TRUE, will undertake statistical filtering based on the
#' a p-value or adjusted p-value threshold stated by `alpha`.Default
#' set at FALSE. Requires presence of columns containing statistical data.
#' In order to filter by the adjusted p-value, a column named `padjusted` must
#' be present. See [mobileRNA::RNAdifferentialAnalysis()] to calculate
#' statistical values.
#'
#' @param alpha numeric; adjusted p-value cutoff as the target FDR for
#' independent filtering. Default is 0.1. Only mobile molecules with adjusted
#' p-values equal or lower than specified are returned.
#'
#'@param threshold numeric; set a threshold level. For sRNA analysis, this
#'represents filtering by the minimum number of replicates that defined the
#' consensus dicercall which is stored in the `DicerCounts` column. While,
#' for mRNA analysis this represents the number of replicates which contained
#' reads for the mRNA molecule which is stored in the `SampleCounts` column.
#'
#'@param baymobil logical; whether to utilise baymobil computation. Only
#'suitable for mRNA. Default is FALSE.
#'
#'@param directory path; path to stored .bam alignment files. For baymobil only.
#'
#'@param heterografts character vector; containing names of heterograft samples.
#'as shown on .bam alignment files.For baymobil only.
#'
#'@param genomefile path; path to a FASTA genome reference file. For baymobil
#'only.
#'@param output_dir path; directory to store output. For baymobil only.
#'
#'@param nmax numeric; nmax parameter for baymobile. Default is 10. For
#'baymobil only.
#'
#'@param Bayes_cutoff numeric; log10 Bayes factor threshold for mobility. For
#'baymobil only. Default is 1.
#'
#'@param condaenv character; name or directory of the Conda environment to use
#' where OS dependencies are stored. For baymobil only.
#'
#'@param snpEff_database character; name of snpEff database to utilise. If
#'custom snpEff database is utilised, provide given database name. For baymobil
#'only. Default is "Arabidopsis_thaliana".
#'
#'@param custom_database logical; state whether using a custom snpEff database.
#'For baymobil only.
#'@param GenomeA.ID character; string to represent prefix added to
#'existing chromosome names in `genomeA`. Default set as "A". For baymobil only.
#'
#'@param GenomeB.ID character; string to represent prefix added to
#'existing chromosome names in `genomeB`. Default set as "B". For baymobil only.
#'
#'
#' @return A data frame containing candidate mobile sRNAs or mRNAs, which could
#' have been further filtered based on statistical significance and the ability
#' to by-pass the thresholds which determine the number of replicates that
#' defined the consensus dicercall (sRNA) or contributed to reads counts (mRNA).
#'
#' More information on `baymobil` output:
#' When employing this computation, the alignment files must be appropriately
#' pre-processed and formatted accordingly. Hence, in your given output
#' directory three folders are output:
#' - 1_variantcalling
#' - 2_snp_df
#' - 3_baymobil
#'
#' The `1_variantcalling` folder stores information generated during variant
#' calling and annotation. These files are important for analysis.
#'
#' The `2_snp_df` folder stores the files utilised to put into baymobil.
#' Specifically, the individual Dataframe representing the total information
#' from homografts (Nh1_Nh2_snp_counts.csv) and the individual Dataframe
#' representing the total information from heterografts (Nn_snp_counts.csv). The
#' file which contains all information is "Nn_Nh_nh_snp_counts.csv".
#'
#' The `3_baymobil` folder stores the output from baymobil computation. There
#' should be one single file named "results_baymobil.csv".
#'
#'Note that this analysis can only handle a maximum of four replicates per
#'condition.
#'
#' @examples
#'data("sRNA_data")
#'
#'
#' # vector of control names
#' controls <- c("selfgraft_1", "selfgraft_2" , "selfgraft_3")
#'
#' # Locate potentially mobile sRNA clusters associated to tomato, no
#' # statistical analysis
#' mobile_df1 <- RNAmobile(input = "sRNA", data =  sRNA_data,
#' controls = controls, genome.ID = "B_", task = "keep", statistical = FALSE)
#'
#'
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyselect starts_with
#' @importFrom dplyr case_when
RNAmobile <- function(input = c("sRNA", "mRNA"), data,
                      controls, genome.ID,
                      task = NULL,
                      statistical = FALSE,
                      alpha = 0.1,
                      threshold = NULL,
                      baymobil = FALSE,
                      directory = NULL,
                      heterografts = NULL,
                      genomefile = NULL,
                      output_dir = NULL,
                      nmax =10,
                      Bayes_cutoff = 1,
                      condaenv = NULL,
                      snpEff_database = "Arabidopsis_thaliana",
                      custom_database = FALSE,
                      GenomeA.ID = "A",
                      GenomeB.ID = "B"){
  if (base::missing(input)) {
    stop("Please specify a character vector of either `sRNA` or `mRNA`
                 to input parameter.")
  }

  if (!base::inherits(data, c("matrix", "data.frame", "DataFrame"))) {
    stop("data must be an object of class matrix, data.frame, DataFrame")
  }
  if (base::missing(controls) || !base::inherits(controls, "character")) {
    stop("Please specify a character vector storing names of control replicates")
  }
  if (base::missing(genome.ID) || genome.ID %in% "") {
    stop("Please specify a single character string which is present in
          the all the chromosomes within the genome you wish to keep or remove")
  }

  if (baymobil == FALSE) {
    y <- data %>%
      dplyr::filter(dplyr::case_when(
        is.null(task) & base::grepl(genome.ID, chr) ~ TRUE,
        task == "remove" & !base::grepl(genome.ID, chr) ~ TRUE,
        task == "keep" & base::grepl(genome.ID, chr) ~ TRUE,
        TRUE ~ FALSE
      ))
    res <- .remove_mapping_errors(data = y, controls = controls)
    if (statistical) {
      res <- res %>% dplyr::filter(padjusted <= alpha)
    }
    if(!is.null(threshold)){
      if(input == "sRNA"){
        res <- res %>% filter(!DicerCounts < threshold)
      }
      if(input == "mRNA"){
        res <- res %>% filter(!SampleCounts < threshold)
      }
    }

    # remove zero values
    zero_count_rows <- rowSums(res[grep("^Count_", names(res))] == 0) == sum(grepl("^Count_", names(res)))

    # Subset the dataframe to remove rows with all zero values in "Count_" columns
    res_fin <- res[!zero_count_rows, ]

    return(res_fin)
  } else
    if(baymobil == TRUE){

      # connect to conda
      reticulate::use_condaenv(condaenv, required = TRUE)
      t <- reticulate::py_config()
      exists_res <- t$prefix
      if(exists_res != condaenv){
        stop("R has not connected to the conda environment")
      }

      # check if software is installed bcftools + snpEFF + baymobil

      t <- snpEFF_exists()
      if(length(t) > 0) {
        message("snpEff package located in conda env...")
      } else {
          stop("Could not locate snpEff in conda env.") }

      b <- bcftools_exists()
      if(length(b) > 0) {
        message("bcftools package located in conda env...")
      } else {
        stop("Could not locate bcftools in conda env.") }


      # set vars
      var1 <- genomefile
      var2 <- paste0(directory,"/", controls, ".bam")
      var3 <- paste0(directory,"/", heterografts, ".bam")
      var4 <- output_dir
      var5 <- nmax
      var6 <- snpEff_database
      var7 <- GenomeA.ID
      var8 <- GenomeB.ID
      var9 <- custom_database

      # check for algnment files.
      bamFiles <- file.exists(c(var2, var3))
      if (bamFiles == FALSE) {
        stop("ERROR: alignment files not found.")
      }

      # check fasta exists
      fasta_exist <- file.exists(genomefile)
      if (fasta_exist == FALSE) {
        stop("ERROR: genome fasta file not found in directory.")
      }

      # Path to the bash script
      baymobil_path <- system.file("baymobil","baymobil_run.sh",
                                   package="mobileRNA")

      # Construct the command
      command <- sprintf("bash %s %s %s &",
                         baymobil_path, var1, var2, var3,var4, var5, var6,
                         var7, var8, var9 )

      message("Initiating analysis with baymobil...")
      message("Initiating analysis with baymobil...")

      # Execute the command in the background
      system(command, wait = FALSE)
      message("Analysis with baymobil complete. Completing processing")

      #import results
      res <- read.csv(paste0(output_dir, "/4_baymobil/results_baymobil.csv"))
      annotaton <- read.csv(paste0(output_dir,"output_dir/1_variantcalling/4_Heterografts_variants_reformatted.txt"))
      combined <- merge(res, annotaton, by = "SNP")

      # select obly positive values
      mobile <- res[res$log10BF > Bayes_cutoff,]

      # rank & label
      mobile$BF <-paste0(round(10^(mobile$log10BF)), "x")
      mobile$confidence_rank <- rank(-mobile$log10BF )

      # filter original df, and bind results.

      # return output
      return()
    }


}
