## code to prepare `res_sRNA_24_DESeq` dataset goes here

groups <- c("Tomato/Eggplant", "Tomato/Eggplant", "Tomato/Eggplant",
            "Tomato/Tomato", "Tomato/Tomato", "Tomato/Tomato")


sRNA_24_edgeR <- RNAanalysis(data = sRNA_24,
                             group = groups,
                             method = "edgeR" )

usethis::use_data(sRNA_24_edgeR, overwrite = TRUE)

