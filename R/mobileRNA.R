#' mobileRNA: Explore RNA mobilome & population-scale changes
#'
#' Uses small RNA or messenger RNA sequencing data in two conditions and 
#' identifies changes in the RNA population, but mainly, identify a putative RNA 
#' mobilome in a chimeric system. For example, in plant graft systems. The 
#' input is the report files produced by the pre-processing steps and the output 
#' is a data-frame of RNA molecules. Output can also consist of the consensus 
#' RNA sequence of small RNA clusters within the determined population. Includes
#' quality control and plotting functions. 
#'
#' The most important functions in the \pkg{mobileRNA} are:
#' \describe{
#' \item{\code{\link{RNAmergeGenomes}}}{Merge two genome assembly files (FASTA)}
#'  \item{\code{\link{RNAmergeAnnotations}}}{Merge two genome annotation files
#'  (GFF)}
#'  \item{\code{\link{mapRNA}}}{Pre-processing of sRNAseq and mRNAseq 
#'  (alignment, raw count, cluster analysis)}
#'  \item{\code{\link{RNAimport}}}{Reads the pre-processing report files in 
#'  to a dataframe for all conditions}
#'  \item{\code{\link{RNAdicercall}}}{Calculates the consensus sRNA 
#'  dicercall class}
#'  \item{\code{\link{RNAsubset}}}{Subsets the working data set to select a 
#'  defined sRNA class}
#'  \item{\code{\link{RNAdifferentialAnalysis}}}{Undertakes differential 
#'  analysis using either the edgeR or DESeq2 method for data structure}
#'  \item{\code{\link{RNAmobile}}}{Identifies mobile candidates in chimeric 
#'  systems}
#'  \item{\code{\link{RNApopulation}}}{Selects unique RNAs within the defined
#'  sample replicates}
#'  \item{\code{\link{RNAsummary}}}{Summarise the differential abundance of RNAs}
#'   \item{\code{\link{RNAreorder}}}{Reorder the data frame for differential 
#'   analysis, ensuring control verse treatment comparison}
#'  \item{\code{\link{RNAsequences}}}{Defines the consensus RNA sequence for a 
#'  sRNA cluster}
#'  \item{\code{\link{RNAattributes}}}{Overlaps the sRNA clusters with a GFF
#'  annotation file and adds overlapping features to the sRNA cluster 
#'  information}
#'  \item{\code{\link{RNAdistribution}}}{Plot the distribution of sRNA 
#'  lengths/classes}
#'  \item{\code{\link{plotHeatmap}}}{Plots a heatmap of log2-transformed 
#'  normalised RPM/FPKM values}
#'  \item{\code{\link{plotSampleDistance}}}{Plots a sample distance heatmap for
#'  quality control}
#'  \item{\code{\link{plotSamplePCA}}}{Plots a PCA plot, customise ratio, 
#'  colours and shapes}
#'  \item{\code{\link{RNAfeatures}}}{Overlaps the sRNA clusters with a GFF
#'  annotation file and calculates the percentage of clusters which overlap with 
#'  specific genomic features}
#'  \item{\code{\link{RNAmean}}}{Calculate mean RPM and mean raw count across 
#'  given samples}
#'  \item{\code{\link{RNAdf2se}}}{Convert `mobileRNA` dataframe to 
#'  a SummarizedExperiment object }
#' }
#'
#' @author
#' Katie Jeynes-Cupper \email{kej031@@student.bham.ac.uk},
#' Marco Catoni \email{m.catoni@@bham.ac.uk}
#' Maintainer: Katie Jeynes-Cupper \email{kej031@@student.bham.ac.uk}
#' @name mobileRNA
#' @docType package
#' @seealso See \code{vignette("mobileRNA", package = "mobieRNA")} for an 
#' overview of the package.
NULL


#' @name sRNA_data
#' @title sRNA_data: simulated data for biological replicates
#' @docType data
#' @usage data(sRNA_data)
#' @description Simulated sRNAseq dataset
#' @return Dataframe in global environment 
#' @details Simulates data is taken from eggplant and tomato sRNAseq samples and
#' created to simulate to movement of sRNA molecules from an Tomato rootstock to
#' an Eggplant Scion. Three Eggplant replicates were spiked with the same 150
#' tomato sRNA clusters, and named "heterograft_" 1 to 3. The analysis compares
#' these heterografts to three Eggplant self-grafts which are the original
#'  un-spiked Eggplant replicates, called "selfgraft_" 1 to 3. 
#'  
#'  This data was imported and organised by the [RNAimport()] function. 
#' @examples
#'  data("sRNA_data")
NULL


#' @name mRNA_data
#' @title mRNA_data: simulated messenger RNA data for biological replicates
#' @docType data
#' @usage data(mRNA_data)
#' @description Simulated mRNAseq dataset
#' @return Dataframe in global environment 
#' @details Simulates data is taken from eggplant and tomato mRNAseq samples and
#' created to simulate to movement of mRNA molecules from an Tomato rootstock to
#' an Eggplant Scion. Two Eggplant replicates were spiked with the same 150
#' tomato mRNA clusters, and named "heterograft_" 1 to 2. The analysis compares
#' these heterografts to two Eggplant self-grafts which are the original
#'  un-spiked Eggplant replicates, called "selfgraft_" 1 to 2. 
#'  
#'  This data was imported and organised by the [RNAimport()] function. 
#' @examples
#'  data("mRNA_data")
NULL

#' @name sRNA_data_dicercall
#' @title sRNA_data_dicercall: simulated data for biological replicates
#' @docType data
#' @keywords dataset
#' @usage data(sRNA_data_dicercall)
#' @description Simulated sRNAseq dataset with RNA consensus
#' @details Simulates data is taken from eggplant and tomato sRNAseq samples and
#' created to simulate to movement of sRNA molecules from an Tomato rootstock to
#' an Eggplant Scion. Three Eggplant replicates were spiked with the same 150
#' tomato sRNA clusters, and named "heterograft_" 1 to 3. The analysis compares
#' these heterografts to three Eggplant self-grafts which are the original
#' un-spiked Eggplant replicates,  called "selfgraft_" 1 to 3.
#' 
#' This data is a further development of the `sRNA_data` data set, where the 
#' [RNAdicercall()] function has been employed. 
#' @return Dataframe in global environment 
#' @examples
#'  data("sRNA_data_dicercall")
NULL

#' @name sRNA_data_mobile
#' @title sRNA_data_mobile: simulated data for biological replicates
#' @docType data
#' @keywords dataset
#' @usage data(sRNA_data_mobile)
#' @description Simulated sRNAseq dataset - potentially mobile RNAs
#' @return Dataframe in global environment 
#' @details Simulates data is taken from eggplant and tomato sRNAseq samples and
#' created to simulate to movement of sRNA molecules from an Tomato rootstock to
#' an Eggplant Scion. Three Eggplant replicates were spiked with the same 150
#' tomato sRNA clusters, and named "heterograft_" 1 to 3. The analysis compares
#' these heterografts to three Eggplant self-grafts which are the original
#' un-spiked Eggplant replicates, called "selfgraft_" 1 to 3. 
#' 
#' This data is a further development of the `sRNA_data` data set, where the
#' candidate mobile sRNAs have been selected. 
#' @examples
#'  data("sRNA_data_mobile")
NULL



#' @name reduced_chr12_Eggplant.gff
#' @title reduced_chr12_Eggplant.gff
#' @docType data
#' @keywords GFF
#' @usage base::system.file("extdata","reduced_chr12_Eggplant.gff",package="mobileRNA")
#' @description Sample GFF file 
#' @return GFF file 
#' @details A small selection of genomic features from chromosome 12 of 
#' Eggplant (Solanum melongena) from Sol Genomics Network 
#' (https://solgenomics.net/)
#' @examples
#'  reduced_chr12_Eggplant <- base::system.file("extdata",
#'  "reduced_chr12_Eggplant.gff",package="mobileRNA")
NULL


#' @name reduced_chr2_Tomato.gff
#' @title reduced_chr2_Tomato.gff
#' @docType data
#' @keywords GFF
#' @usage base::system.file("extdata","reduced_chr2_Tomato.gff",package="mobileRNA")
#' @description Sample GFF file 
#' @return GFF file 
#' @details A small selection of genomic features from chromosome 2 of 
#' Tomato (Solanum lycopersicum) from Sol Genomics Network 
#' (https://solgenomics.net/)
#' @examples
#'  reduced_chr2_Tomato <- base::system.file("extdata",
#'  "reduced_chr2_Tomato.gff",package="mobileRNA")
NULL


#' @name reduced_chr12_Eggplant.fa
#' @title reduced_chr12_Eggplant.fa
#' @docType data
#' @keywords FASTA
#' @usage base::system.file("extdata","reduced_chr12_Eggplant.fa",package="mobileRNA")
#' @description Sample FASTA file 
#' @return FASTA file 
#' @details Reduced nucleotide sequence of chromosome 12 of 
#' Eggplant (Solanum melongena) from Sol Genomics Network 
#' (https://solgenomics.net/)
#' @examples
#'  reduced_chr12_Eggplant <- base::system.file("extdata",
#'  "reduced_chr12_Eggplant.fa",package="mobileRNA")
NULL


#' @name reduced_chr2_Tomato.fa
#' @title reduced_chr2_Tomato.fa
#' @docType data
#' @keywords FASTA
#' @usage base::system.file("extdata","reduced_chr2_Tomato.fa",package="mobileRNA")
#' @description Sample FASTA file 
#' @return FASTA file 
#' @details Reduced nucleotide sequence of chromosome 2 of Tomato 
#' (Solanum lycopersicum) from Sol Genomics Network 
#' (https://solgenomics.net/)
#' @examples
#'  reduced_chr2_Tomato <- base::system.file("extdata",
#'  "reduced_chr2_Tomato.fa",package="mobileRNA")
NULL


#' @name selfgraft_demo_1.fq
#' @title selfgraft_demo_1.fq
#' @docType data
#' @keywords FASTQ
#' @usage base::system.file("extdata/sRNAseq","selfgraft_demo_1.fq",package="mobileRNA")
#' @description Demo FASTQ file 
#' @return FASTQ file 
#' @details Demo sRNAseq FASTQ file for a self-graft control Eggplant 
#' (Solanum melongena) sample. This file has been majorly reduced and does not
#' represent true data. For pre-processing purposes. 
#' @examples
#'  selfgraft_demo_1 <- base::system.file("extdata/sRNAseq","selfgraft_demo_1.fq",
#'  package="mobileRNA")
NULL

#' @name selfgraft_demo_2.fq
#' @title selfgraft_demo_2.fq
#' @docType data
#' @keywords FASTQ
#' @usage base::system.file("extdata/sRNAseq","selfgraft_demo_2.fq",package="mobileRNA")
#' @description Demo FASTQ file 
#' @return FASTQ file 
#' @details Demo sRNAseq FASTQ file for a self-graft control Eggplant 
#' (Solanum melongena) sample. This file has been majorly reduced and does not
#' represent true data. For pre-processing purposes. 
#' @examples
#'  selfgraft_demo_2 <- base::system.file("extdata/sRNAseq","selfgraft_demo_2.fq",
#'  package="mobileRNA")
NULL

#' @name selfgraft_demo_3.fq
#' @title selfgraft_demo_3.fq
#' @docType data
#' @keywords FASTQ
#' @usage base::system.file("extdata/sRNAseq","selfgraft_demo_3.fq",package="mobileRNA")
#' @description Demo FASTQ file 
#' @return FASTQ file 
#' @details Demo sRNAseq FASTQ file for a self-graft control Eggplant 
#' (Solanum melongena) sample. This file has been majorly reduced and does not
#' represent true data. For pre-processing purposes. 
#' @examples
#'  selfgraft_demo_3 <- base::system.file("extdata/sRNAseq","selfgraft_demo_3.fq",
#'  package="mobileRNA")
NULL


#' @name heterograft_demo_1.fq
#' @title heterograft_demo_1.fq
#' @docType data
#' @keywords FASTQ
#' @usage base::system.file("extdata/sRNAseq","heterograft_demo_1.fq",package="mobileRNA")
#' @description Demo FASTQ file 
#' @return FASTQ file 
#' @details Demo sRNAseq FASTQ file for a heterograft sample where the sample 
#' was taken from the shoot of the scion. The scion is Eggplant 
#' (Solanum melongena) and the rootstock is Tomato (Solanum lycopersicum). 
#' This file has been majorly reduced and does not represent true data. 
#' For pre-processing purposes. 
#' @examples
#'  heterograft_demo_1 <- base::system.file("extdata/sRNAseq","heterograft_demo_1.fq",
#'  package="mobileRNA")
NULL

#' @name heterograft_demo_2.fq
#' @title heterograft_demo_2.fq
#' @docType data
#' @keywords FASTQ
#' @usage base::system.file("extdata/sRNAseq","heterograft_demo_2.fq",package="mobileRNA")
#' @description Demo FASTQ file 
#' @return FASTQ file 
#' @details Demo sRNAseq FASTQ file for a heterograft sample where the sample 
#' was taken from the shoot of the scion. The scion is Eggplant 
#' (Solanum melongena) and the rootstock is Tomato (Solanum lycopersicum). 
#' This file has been majorly reduced and does not represent true data. 
#' For pre-processing purposes. 
#' @examples
#'  heterograft_demo_2 <- base::system.file("extdata/sRNAseq","heterograft_demo_2.fq",
#'  package="mobileRNA")
NULL

#' @name heterograft_demo_3.fq
#' @title heterograft_demo_3.fq
#' @docType data
#' @keywords FASTQ
#' @usage base::system.file("extdata/sRNAseq","heterograft_demo_3.fq",package="mobileRNA")
#' @description Demo FASTQ file 
#' @return FASTQ file 
#' @details Demo sRNAseq FASTQ file for a heterograft sample where the sample 
#' was taken from the shoot of the scion. The scion is Eggplant 
#' (Solanum melongena) and the rootstock is Tomato (Solanum lycopersicum). 
#' This file has been majorly reduced and does not represent true data. 
#' For pre-processing purposes. 
#' @examples
#'  heterograft_demo_3 <- base::system.file("extdata/sRNAseq","heterograft_demo_3.fq",
#'  package="mobileRNA")
NULL

#' @name selfgraft_mRNAdemo_1.fq
#' @title selfgraft_mRNAdemo_1.fq
#' @docType data
#' @keywords FASTQ
#' @usage base::system.file("extdata/mRNAseq","selfgraft_mRNAdemo_1.fq",package="mobileRNA")
#' @description Demo FASTQ file 
#' @return FASTQ file 
#' @details Demo mRNAseq FASTQ file for a self-graft control Eggplant 
#' (Solanum melongena) sample. This file has been majorly reduced and does not
#' represent true data. For pre-processing purposes. 
#' @examples
#'  selfgraft_mRNAdemo_1 <- base::system.file("extdata/mRNAseq","selfgraft_mRNAdemo_1.fq",
#'  package="mobileRNA")
NULL

#' @name selfgraft_mRNAdemo_2.fq
#' @title selfgraft_mRNAdemo_2.fq
#' @docType data
#' @keywords FASTQ
#' @usage base::system.file("extdata/mRNAseq","selfgraft_mRNAdemo_2.fq",package="mobileRNA")
#' @description Demo FASTQ file 
#' @return FASTQ file 
#' @details Demo mRNAseq FASTQ file for a self-graft control Eggplant 
#' (Solanum melongena) sample. This file has been majorly reduced and does not
#' represent true data. For pre-processing purposes. 
#' @examples
#'  selfgraft_mRNAdemo_2 <- base::system.file("extdata/mRNAseq","selfgraft_mRNAdemo_2.fq",
#'  package="mobileRNA")
NULL


#' @name heterograft_mRNAdemo_1.fq
#' @title heterograft_mRNAdemo_1.fq
#' @docType data
#' @keywords FASTQ
#' @usage base::system.file("extdata/mRNAseq","heterograft_mRNAdemo_1.fq",package="mobileRNA")
#' @description Demo FASTQ file 
#' @return FASTQ file 
#' @details Demo mRNAseq FASTQ file for a heterograft sample where the sample 
#' was taken from the shoot of the scion. The scion is Eggplant 
#' (Solanum melongena) and the rootstock is Tomato (Solanum lycopersicum). 
#' This file has been majorly reduced and does not represent true data. 
#' For pre-processing purposes. 
#' @examples
#'  heterograft_mRNAdemo_1 <- base::system.file("extdata/mRNAseq",
#'  "heterograft_mRNAdemo_1.fq", package="mobileRNA")
NULL


#' @name heterograft_mRNAdemo_2.fq
#' @title heterograft_mRNAdemo_2.fq
#' @docType data
#' @keywords FASTQ
#' @usage base::system.file("extdata/mRNAseq","heterograft_mRNAdemo_2.fq",package="mobileRNA")
#' @description Demo FASTQ file 
#' @return FASTQ file 
#' @details Demo mRNAseq FASTQ file for a heterograft sample where the sample 
#' was taken from the shoot of the scion. The scion is Eggplant 
#' (Solanum melongena) and the rootstock is Tomato (Solanum lycopersicum). 
#' This file has been majorly reduced and does not represent true data. 
#' For pre-processing purposes. 
#' @examples
#'  heterograft_mRNAdemo_2 <- base::system.file("extdata/mRNAseq",
#'  "heterograft_mRNAdemo_2.fq", package="mobileRNA")
NULL



