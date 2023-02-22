## code to prepare `res_sRNA_24_DESeq` dataset goes here

groups <- c("Tomato/Eggplant", "Tomato/Eggplant", "Tomato/Eggplant",
            "Tomato/Tomato", "Tomato/Tomato", "Tomato/Tomato")


## Differential analysis: DEseq2 method
# 24-nt sRNA data-set
sRNA_24_DESeq2 <- RNAanalysis(data = sRNA_24,
                              group = groups,
                              method = "DESeq2" )


usethis::use_data(sRNA_24_DESeq2, overwrite = TRUE)






