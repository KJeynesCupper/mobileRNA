########### ########### ########### ###########
########### Adding the internal data ###########
########### ########### ########### ###########


########### ########### ########### ########### ########### ########### ###########


## SUMMARY DATA.
## Path to genomic overlap data file (.bed)
genomic_features <-  "/Users/U2093048/Projects/relocate/build_internal_data/locifile_annotation_overlap_TomEgg.bed"

## Import & organise data.
results_dir <-  "/Users/U2093048/Projects/relocate/build_internal_data/second_alignment_unique_mincov/"
conditions <- c("TomEgg_1","TomEgg_2", "TomEgg_3",
                "TomTom_1", "TomTom_2", "TomTom_3")
totalNumReads <- c(24441759,21378845, 21482356, 3951725, 3343954, 2586910)
names(totalNumReads) <- conditions
clusterlocations <- rtracklayer::import.gff("/Users/U2093048/Projects/relocate/build_internal_data/gff_loci_TomEgg_vs_Tom.gff3")

sRNA_data <- sample_table(results_dir,
                           conditions,
                           clusterlocations,
                           totalNumReads,
                           genomic_features)

#### Save data
usethis::use_data(sRNA_data)



########### ########### ########### ########### ########### ########### ###########

##  define consensus sRNA classes.
samples <- c("TomEgg_1", "TomEgg_2", "TomEgg_3")

sRNA_data_summary <- define_consensus(data = sRNA_data,
                                 samples = samples,
                                 tidy=TRUE)

# Subset data: 24-nt siRNAs
sRNA_24 <- class_subset(sRNA_data_summary, type = 24)
# Subset data: 24 21/22-nt siRNAs
sRNA_2122 <- class_subset(sRNA_data_summary, type = c(21, 22))


########### ########### ########### ########### ########### ########### ###########

## DESeq Normalisation for each subset.

# sample conditions.
group <- c("Tomato/Eggplant", "Tomato/Eggplant", "Tomato/Eggplant",
           "Tomato/Tomato", "Tomato/Tomato", "Tomato/Tomato")

# 24-nt subset
sRNA_24_norm_DESeq <- normalise_data(sRNA_24, group, method = "DESeq2" )

# 21/22-nt subset
sRNA_2122_norm_DESeq <- normalise_data(sRNA_2122, group, method = "DESeq2" )


## edgeR Normalisation for each subset.
# 24-nt subset
sRNA_24_norm_edgeR <- normalise_data(sRNA_24, group, method = "edgeR" )
# 21/22-nt subset
sRNA_2122_norm_edgeR <- normalise_data(sRNA_2122, group, method = "edgeR" )



########### ########### ########### ########### ########### ########### ###########

## Differential analysis: DEseq2 method
sRNA_24_DE_DESeq <- DE_analysis(data = sRNA_24_norm_DESeq,
                                groupPair = c("Tomato/Eggplant","Tomato/Tomato"),
                                method = "DESeq2" )

sRNA_2122_DE_DESeq <- DE_analysis(data = sRNA_2122_norm_DESeq,
                                  groupPair = c("Tomato/Eggplant","Tomato/Tomato"),
                                  method = "DESeq2" )

## Add results to the summary data frame.
res_sRNA_24_DESeq <- add_results(sRNA_24,sRNA_24_DE_DESeq, method = "DESeq2")
res_sRNA_2122_DESeq <- add_results(sRNA_2122,sRNA_2122_DE_DESeq, method = "DESeq2")


## Differential analysis: edgeR method
sRNA_24_DE_edgeR <- DE_analysis(data = sRNA_24_norm_edgeR,
                                groupPair = c("Tomato/Eggplant","Tomato/Tomato"),
                                method = "edgeR" )

sRNA_2122_DE_edgeR <- DE_analysis(data = sRNA_2122_norm_edgeR,
                                  groupPair = c("Tomato/Eggplant","Tomato/Tomato"),
                                  method = "edgeR" )

## Add results to the summary data frame.
res_sRNA_24_edgeR <- add_results(sRNA_24,sRNA_24_DE_edgeR, method = "edgeR")
res_sRNA_2122_edgeR <- add_results(sRNA_2122,sRNA_2122_DE_edgeR, method = "edgeR")


########### ########### ########### ########### ########### ########### ###########



## identify mobile molecules


# vector of control names
controls <- c("TomTom_1", "TomTom_2", "TomTom_3")

## Mobile Molecules: DEseq2 method

# remove clusters associated to tomato in 24nt RNA database
sRNA_24_mobile_DEseq <- mobile_molecules(data = res_sRNA_24_DESeq, controls, id = "SL40", task = "remove")

# remove clusters associated to tomato in 2122nt RNA database
sRNA_2122_mobile_DEseq  <- mobile_molecules(data = res_sRNA_2122_DESeq,controls, id = "SL40", task = "remove")



## Mobile Molecules: edgeR method

# remove clusters associated to tomato in 24nt RNA database
sRNA_24_mobile_edgeR <- mobile_molecules(data = res_sRNA_24_edgeR, controls, id = "SL40", task = "remove")

# remove clusters associated to tomato in 2122nt RNA database
sRNA_2122_mobile_edgeR <- mobile_molecules(data = res_sRNA_2122_edgeR,controls, id = "SL40", task = "remove")


########### ########### ########### ########### ########### ########### ###########



