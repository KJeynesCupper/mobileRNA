## code to prepare `sRNA_data_summary` dataset goes here
samples <- c("TomEgg_1", "TomEgg_2", "TomEgg_3")

# call define_consensus from package
sRNA_data_summary <- RNAconsensus(data = sRNA_data,
                                      conditions = samples,
                                      tidy=TRUE)

usethis::use_data(sRNA_data_summary, overwrite = TRUE)





##### new 06/05/2023
sRNA_data_consensus <- mobileRNA::RNAconsensus(sRNA_data,
                                    conditions = c("heterograft_1", "heterograft_2", "heterograft_3"))

usethis::use_data(sRNA_data_consensus, overwrite = TRUE)

##### new 06/05/2023
sRNA_data_mobile <- mobileRNA::RNAmobile(data = sRNA_data_consensus, controls = c("selfgraft_1", "selfgraft_2" , "selfgraft_3"), id = "SL", task ="keep", statistical = FALSE )
usethis::use_data(sRNA_data_mobile, overwrite = TRUE)
