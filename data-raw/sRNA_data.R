## code to prepare `sRNA_data` dataset goes here
## Path to genomic overlap data file (.bed)
#genomic_features <-  "/Users/U2093048/Projects/package-backup/edits/build_internal_data/locifile_annotation_overlap_TomEgg.bed"

## Import & organise data.
results_dir <-  "../build-data/"

# sample names
conditions <- c("TomEgg_1","TomEgg_2", "TomEgg_3",
                "TomTom_1", "TomTom_2", "TomTom_3")


# import gff file of sRNA cluster info
clusterlocations <- rtracklayer::import.gff("../gff_loci_TomEgg_vs_Tom.gff3")

# call sample_table function to make data.
sRNA_data <- RNAimport(results = results_dir,
                          samples = conditions,
                          clusters = clusterlocations)

#### Save data
usethis::use_data(sRNA_data, overwrite = TRUE)
###############

#updated version of import : synthetic data - 06/07/23
data <- read.delim("../sRNA_synthetic_catoni_data.txt", header = TRUE, sep = "\t", dec = ".")

df2 <- data[order(data$chr),]

Cluster <-  paste0("cluster_", 1:nrow(df2))
sRNA_data <- as.data.frame(append(df2, list(Cluster = Cluster), after = 4))


# name columnns
library(stringr)
library(dplyr)

sRNA_data <- sRNA_data %>%
  rename_with(~str_replace(., 'spiked_eggplant_1', 'heterograft_1'))%>%
  rename_with(~str_replace(., 'spiked_eggplant_2', 'heterograft_2'))%>%
  rename_with(~str_replace(., 'spiked_eggplant_3', 'heterograft_3'))%>%
  rename_with(~str_replace(., 'eggplant_1', 'selfgraft_1'))%>%
  rename_with(~str_replace(., 'eggplant_2', 'selfgraft_2'))%>%
  rename_with(~str_replace(., 'eggplant_3', 'selfgraft_3'))

usethis::use_data(sRNA_data, overwrite = TRUE, compress = "xz")

tools::resaveRdaFiles("./data/",compress="xz")
tools::checkRdaFiles("data/")
