## code to prepare `sRNA_2122_mobile_DEseq` dataset goes here

## identify mobile molecules

# vector of control names
controls <- c("TomTom_1", "TomTom_2", "TomTom_3")

## Mobile Molecules: DEseq2 method

# remove clusters associated to tomato in 2122nt RNA database
sRNA_2122_mobile_DESeq2  <- RNAmobile(data = sRNA_2122_DESeq2,controls, id = "SL40", task = "remove")



usethis::use_data(sRNA_2122_mobile_DESeq2, overwrite = TRUE)
