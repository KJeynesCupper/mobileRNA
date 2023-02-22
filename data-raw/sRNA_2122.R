## code to prepare `sRNA_2122` dataset goes here

# this creates a dataframe from summary data which only contains sRNA with a length og either 21 or 22 nucleotides.
sRNA_2122 <- RNAsubset(sRNA_data_summary, type = c(21, 22), sig = FALSE)

usethis::use_data(sRNA_2122, overwrite = TRUE)
