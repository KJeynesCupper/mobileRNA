#------------------------------------------------------------#
# Title:  Invisible functions                                #
# Author: Katie Jeynes-Cupper (katie.jeynescupper@gmail.com) #
# Date:   01.02.23                                           #
#------------------------------------------------------------#



################ remove mapping errors (RNAmobile function) ####################
.remove_mapping_errors <- function(data, controls) {
  class_colnames  <- data %>% dplyr::select(paste0("Count_", controls))

  if (length(colnames(class_colnames)) > 1){
    x <- c()
    for (j in 1:nrow(data)){
      if(sum(stats::na.omit(as.numeric( data[j,colnames(class_colnames)],
                                        na.rm=TRUE)))>0){
        x <- c(x,j)
      }
    }
  } else
    if (length(colnames(class_colnames)) == 1){
      x <- c()
      for (k in 1:nrow(data)){
        if(stats::na.omit(as.numeric(data[k,colnames(class_colnames)],
                                     na.rm=TRUE))!= 0){
          x <- c(x,k)
        }
      }
    }
  if(is.null(x)){
    data <- data
  } else
  data <- data[-x,]

  return(data)
}


################ DESE2 function (RNAanalysis function) #########################
.DESeq_normalise <- function(data, conditions){
  column.data <- data.frame(conditions=as.factor(conditions))
  base::rownames(column.data) <- base::colnames(data)
  count.data.set <- DESeq2::DESeqDataSetFromMatrix(countData=data,
                                                   colData=column.data,
                                                   design= ~ conditions)
  count.data.set$conditions <- stats::relevel(count.data.set$conditions,
                                              conditions[1])
  out <- DESeq2::estimateSizeFactors(count.data.set)
  return(out)
}

################ EDGER function (RNAanalysis function) #########################
.edgeR_normalise <- function(data, conditions){
  d <- edgeR::DGEList(counts = data, group = factor(conditions))
  result <- edgeR::calcNormFactors(d)
  result$samples
  result<- edgeR::estimateDisp(result)
  result$common.dispersion
  return(result)
}

################ Find RNA complementary sequence (RNAsequence function) ########

find_complementary_sequenceRNA <- function(seq) {
  # conversions
  conversion_nucleotides <- c(A = "U", U = "A", C = "G", G = "C")

  # calculate complementary nt for each
  complementary_calc <- sapply(strsplit(seq, ""), function(nucleotide) {
    conversion_nucleotides[nucleotide]
  })

  # Combine into string
  output <- paste0(complementary_calc, collapse = "")
  return(output)
}


################ Find DNA complementary sequence (RNAsequence function) #######

find_complementary_sequenceDNA <- function(seq) {
  # conversions
  conversion_nucleotides <- c(A = "T", U = "A", C = "G", G = "C")

  # calc complementary nt for each in string
  complementary_calc <- sapply(strsplit(seq, ""), function(nucleotide) {
    conversion_nucleotides[nucleotide]
  })

  # Combine into string
  output <- paste0(complementary_calc, collapse = "")
  return(output)
}


################## global variable storage #####################################
utils::globalVariables(c("ID", "sRNA_Consensus", "nt_20", "nt_21", "nt_22",
                         "nt_23", "nt_24", "group", "qc", "score", "phase",
                         "type", "V1", "V10", "V11", "V12", "V13", "V14", "V15",
                         "V16", "V17", "V2", "V3", "V4", "V5", "V6", "V7", "V8",
                         "V9", "padj",".", "nt_N", "n", "FPKM_mean", "chr",
                         "key", "Count", "Class", "padjusted", "pvalue", "freq",
                         "value" , "variable" , "repeats_info" , "Genome" ,
                         "Dataset" ,"setNames" , "DicerCall" , "Reads" , "RPM" ,
                         "MajorRNA", "i"))


