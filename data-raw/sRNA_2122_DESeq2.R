## code to prepare `res_sRNA_24_DESeq` dataset goes here

groups <- c("Tomato/Eggplant", "Tomato/Eggplant", "Tomato/Eggplant",
            "Tomato/Tomato", "Tomato/Tomato", "Tomato/Tomato")

# 2122-nt sRNA data-set
sRNA_2122_DESeq2 <- RNAanalysis(data = sRNA_2122,
                                group = groups,
                                method = "DESeq2" )


usethis::use_data(sRNA_2122_DESeq2, overwrite = TRUE)

