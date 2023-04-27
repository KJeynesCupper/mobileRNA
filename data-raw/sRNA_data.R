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
