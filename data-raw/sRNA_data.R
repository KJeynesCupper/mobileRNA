## code to prepare `sRNA_data` dataset goes here
## Path to genomic overlap data file (.bed)
genomic_features <-  "/Users/U2093048/Projects/relocate/edits/build_internal_data/locifile_annotation_overlap_TomEgg.bed"

## Import & organise data.
results_dir <-  "/Users/U2093048/Projects/relocate/edits/build_internal_data/second_alignment_unique_mincov/"

# sample names
conditions <- c("TomEgg_1","TomEgg_2", "TomEgg_3",
                "TomTom_1", "TomTom_2", "TomTom_3")

# total number of reads from mapping files.
totalNumReads <- c(24441759,21378845, 21482356, 3951725, 3343954, 2586910)
# set names to numbers
names(totalNumReads) <- conditions
# import gff file of sRNA cluster info
clusterlocations <- rtracklayer::import.gff("/Users/U2093048/Projects/relocate/edits/build_internal_data/gff_loci_TomEgg_vs_Tom.gff3")

# call sample_table function to make data.
sRNA_data <- sample_table(results_dir,
                          conditions,
                          clusterlocations,
                          totalNumReads,
                          genomic_features)

#### Save data

usethis::use_data(sRNA_data, overwrite = TRUE)
