## code to prepare `sRNA_24` dataset goes here

# This is a dataframe containing a subset of the summary data, only includinh the sRNAs which have a consensus of 24-nt.
# Subset data: 24-nt siRNAs
sRNA_24 <- RNAsubset(sRNA_data_summary, type = 24, sig = FALSE)
usethis::use_data(sRNA_24, overwrite = TRUE)
