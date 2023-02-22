## code to prepare `sRNA_24_mobile_DEseq` dataset goes here

## identify mobile molecules


# vector of control names
controls <- c("TomTom_1", "TomTom_2", "TomTom_3")

## Mobile Molecules: DEseq2 method

# remove clusters associated to tomato in 24nt RNA database
sRNA_24_mobile_DESeq2 <- mobile_molecules(data = res_sRNA_24_DESeq2, controls, id = "SL40", task = "remove")

usethis::use_data(sRNA_24_mobile_DESeq2, overwrite = TRUE)
