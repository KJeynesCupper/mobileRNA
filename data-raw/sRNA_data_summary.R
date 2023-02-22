## code to prepare `sRNA_data_summary` dataset goes here
samples <- c("TomEgg_1", "TomEgg_2", "TomEgg_3")

# call define_consensus from package
sRNA_data_summary <- RNAconsensus(data = sRNA_data,
                                      conditions = samples,
                                      tidy=TRUE)

usethis::use_data(sRNA_data_summary, overwrite = TRUE)
