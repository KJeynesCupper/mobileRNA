## code to prepare `sRNA_24_mobile_edgeR` dataset goes here

## identify mobile molecules


# vector of control names
controls <- c("TomTom_1", "TomTom_2", "TomTom_3")

## Mobile Molecules: edgeR method

# remove clusters associated to tomato in 24nt RNA database
sRNA_24_mobile_edgeR <- RNAmobile(data = sRNA_24_edgeR, controls, id = "SL40", task = "remove")

usethis::use_data(sRNA_24_mobile_edgeR, overwrite = TRUE)
