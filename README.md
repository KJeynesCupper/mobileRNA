mobileRNA <a href="kjeynescupper.github.io/mobileRNA/"><img src="man/figures/logo.png" align="right" height="138" /></a>
======================================================================
<br>

Overview
======================================================================

**mobileRNA** is an `R` package that provides a pipeline for 
pre-processing and analysis of small RNA (sRNA) and messenger RNA (mRNA) 
sequencing data, primarily for the identification of mobile RNAs in plant graft 
systems. As well as typical treatment vs control analysis to identify changes in
sRNA population such as abundance and production. These two analysis
workflows have be separated as mobile RNA analysis and the core RNA
analysis, respectively.

Plant grafting has been used to study RNA movement and the mechanisms
associated with their action. Research has established that sRNA
molecules can travel between the plant roots and shoots and can
introduce changes to gene expression. There are two key classes of
sRNAs, 24-nt small interfering RNAs (siRNAs) which can introduce DNA
methylation via the RNA-directed DNA methylation pathway and the
21/22-nt microRNAs (miRNAs) which are associated to post-transcriptional
gene silencing. Similarly, studies have shown that mRNAs can move across 
distances, and it is thought they may translate into proteins which act as 
transcription factors in the recipient tissues. Changes in these RNA populations 
could instigate or facilitate grafting-induced traits, such as improved plant 
vigour or stress resistance.

When plants are grafted it joins two distinct genotypes to form a
chimeric plant. The sequence variation between the two genomes involved
can be used to discriminate the origin of a sequenced RNA molecule.
Hence, if an RNA molecule sequenced from tissues of the grafted
partners has found matching the genome of the other grafting partner,
this could empirically demonstrate it's movement across the graft
junction.

Most available genomics approaches to implement this analysis are based
on sRNA sequencing, following by alignment on a genotype of reference
and post alignment screening of genetic variants to identify molecules
which have better match for the genotype of the grafted partner. These
methods have many limitations, which might include:

-   High dependency on arbitrary thresholds to determined mapped and
    unmapped RNA sequencing reads.
-   Use of additional statistical test to discriminate genetic
    variations from sequencing noise in post alignment analysis.
-   High false positive rates.

Here, to circumvents such problems we propose a method inspired by the
RNAseq analysis of plant hybrids, including an alignment step performed
simultaneously on both genomes involved. The rational of this approach
considers that alignment tools already implement an algorithm ideal for
identification of the best matches (accordingly to set parameters) in a
given genome reference, but they do not account for potential matches to
DNA sequences which are not provided as reference. Therefore, the two
genomes from all partners involved in the chimeric system are merged in
a single FASTA file and used as reference for the unique alignment with
the `Bowtie` for sRNAseq or `HISAT2` for mRNAseq. This is in a bid to supply the 
algorithm with as much information as possible to make the best possible
predictions and placement of sequencing reads to each genome. These 
pre-processing steps and downstream analysis is all supported by `mobileRNA`.

The downstream analysis allows for differential analysis for the
comparison in abundance, the identification of changes in RNA
production and the identification of putative mobile RNA molecules via
the filtering of RNA molecules which uniquely map to the mobile genome
of to the other grafting partner involved, which would represent the
putative mobile RNA molecules.

**Below is a quick-start guide to the identification of mobile sRNAs in
a chimeric system** <br>

*Look-out for developments that accommodates mRNA movement* <br>

Author
--------
Katie Jeynes-Cupper, University of Birmingham,
[kej031\@student.bham.ac.uk](mailto:kej031@student.bham.ac.uk){.email}
<br>

Table of Contents
-----------------
-   [Installation](#installation)
-   [Getting Help](#Getting-Help)
-   [Quick Start](#Quick-Start)
-   [Output](#Output)
-   [Advanced Analysis](#Advanced-Analysis)

Installation 
------------------------------------------------------------------------
The latest version of the package can be install directly from this
GitHub repo:

``` r
# Github installation 
if (!require("devtools")) install.packages("devtools")
devtools::install_github("KJeynesCupper/mobileRNA", ref = "main")


# Bioconductor installation 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("mobileRNA")

# load into R library
library(mobileRNA)

```
<br>

Getting Help
------------------------------------------------------------------------
For additional information on each function, please read through the
documentation in the `mobileRNA` package by typing the `?` help operator
before any of the function names in the package or by using the `help()`
function.

For an in-depth step-by-step analysis, consider reading the vignette
provided with this package:

``` r
vignette("mobileRNA")
```

<br>

## Quick Start

Below is the summary workflow:

```{r, fig.align="centre", echo=FALSE,out.width="750", out.height="750",fig.cap="The basic analysis pipeline for mobileRNA" }

knitr::include_graphics("../man/figures/mobileRNA_graphic_1.png")
```

<br><br><br>

## 1. Example Data Set

The package includes a simulated data set to replicate the grafting
between eggplant-tomato where eggplant represents the scion and tomato
represents the rootstock. The FASTQ files represent sRNAseq data
extracted from the eggplant leaf tissue. Here we will locate sRNA
produced by tomato roots which have traveled into the eggplant leaves.

## 2. Merging Genome Assemblies {#merging-genome-assemblies}

Merge the FASTA genome assemblies of tomato and eggplant into a single
reference file stored in your desired directory.

```r
# pull genomes 
library(BiocFileCache)
cache_dir <- tools::R_user_dir("mobileRNA", which = "cache")
cache <- BiocFileCache(cache_dir)

# Construct URL to example FASTA files
url_remote <- "https://github.com/KJeynesCupper/assemblies/raw/main/"

fasta_1_url <- file.path(url_remote, "chr12_Eggplant_V4.1.fa.gz")
fasta_2_url <- file.path(url_remote,"chr2_S_lycopersicum_chromosomes.4.00.fa.gz")

# Download example FASTA files and add them to cache
fasta_1 <- bfcrpath(cache, fasta_1_url)
fasta_2 <- bfcrpath(cache, fasta_2_url)

# define temporary output directory - replace with your directory
output_assembly_file <- file.path(tempfile("merged_annotation", 
                                           fileext = ".fa"))

# merge
merged_reference <- RNAmergeGenomes(genomeA = fasta_1,
                                    genomeB = fasta_2,
                                    out_dir = output_assembly_file)
```

<br>

## 3. Alignment

Align sRNA sequencing reads to the merged genome using our unique
alignment pipeline wrapped by the `mapRNA()` function.

``` r
samples <- file.path(system.file("extdata",package="mobileRNA"))

output_location <- tempdir()

mapRNA(input = "sRNA",
       input_files_dir = samples, 
       output_dir = output_location, 
       genomefile = output_assembly_file,
       condaenv = "ShortStack4",
       mmap = "n")

```

<br>

## 4. Import Pre-Processed Data into R

Import the results from the alignment step into R using the
`RNAimport()` function. This requires the directory storing the sample
output folders and the same of the samples to import from the directory.

``` r
# Directory containing results
results_dir <-  file.path(output_location,"2_alignment_results")

# Sample names and total number of reads, in the same order. 
sample_names <- c("selfgraft_1", "selfgraft_2", "selfgraft_3",
                  "heterograft_1", "heterograft_2", "heterograft_3")


sRNA_data <- RNAimport(input = "sRNA", 
                       directory = results_dir,
                       samples = sample_names)
                           
```

or, load the pre-processed data:

```{r Load, message=FALSE}

data("sRNA_data")

```

<br>

## 5. Calculate the Consensus Dicercall

For a given sRNA cluster, each replicate has determined the dicercall,
also known as the sRNA class, based on the length in nucleotides of the
most abundant sRNA. This can be drawn from all samples or named samples.
The output can be used as threshold values for downstream analysis, and
to remove data noise depending on data quality.

``` r

sRNA_data_summary <- RNAdicercall(data = sRNA_data, tidy = TRUE )
```

<br>

## 6. Differential Analysis of sRNA Population

Undertake differential analysis of sRNA within the experimental design
to explore changes in abundance. The function allows for two methods;
`edgeR` or `DESeq`.

```{r}
## sample conditions in order within dataframe
groups <- c("Selfgraft", "Selfgraft", "Selfgraft", 
            "Heterograft", "Heterograft", "Heterograft")


## Differential analysis of whole dataset: DESeq2 method 
sRNA_DESeq2 <- RNAdifferentialAnalysis(data = sRNA_data,
                                       group = groups,
                                       method = "DESeq2")
                              



                              
# save output as txt file 
write.table(sRNA_DESeq2, "./sRNA_DA_output.txt")

```

Summarise results:

```{r}

RNAsummary(sRNA_DESeq2)

```

How about summarizing the sRNA population which are statistically
significant:

```{r}

RNAsummary(sRNA_DESeq2, alpha=0.05)

```

## 7. Identify Potential Mobile sRNA

Select the putative mobile sRNA clusters using `RNAmobile()`. This
requires supplying the function with a unique identifier of the
rootstock genome. The merging step placed the prefix "B" to the tomato
chromosomes.

```{r}
# define control samples
controls <- c("selfgraft_1", "selfgraft_2", "selfgraft_3")

mobile_sRNA <- RNAmobile(input = "sRNA",
                         data = sRNA_DESeq2, 
                         controls = controls,
                         genome.ID = "B", 
                         task = "keep")



# save output as txt file 
write.table(mobile_sRNA, "./sRNA_mobile_output.txt")

```

<br>

Output
------------------------------------------------------------------------
A data frame where rows represent potential mobile sRNA clusters. The
columns include information on the cluster, individual sample
replicates, and more.

#### Information on the cluster:

-   `Locus`: Name of the chromosome or scaffold, start position & end
    position
-   `chr`: Name of the chromosome or scaffold
-   `start` : Start position of the cluster
-   `end` : End position of the cluster

#### Information on each sample replicate:

-   `Cluster`: Cluster Name
-   `Dicercall` : The size of most abundant small RNA size
-   `Count` : Number of reads, default is uniquely aligned (*e.g.* not
    multi-mapping).
-   `MajorRNA` : RNA sequence of the most abundant sRNA in cluster
    within the sample
-   `RPM` : Reads per Million
-   `FPKM` : Fragments Per Kilobase of transcript per Million

#### Other information

-   `DicerConsensus` : Consensus sRNA class
-   `DicerCounts` : Number of replicates which contributed to the
    consensus dicercall sRNA class
-   `CountMean` : Count mean (Calculated by `RNAdifferentialAnalysis()`)
-   `log2FoldChange` : Log2FoldChange--The effect size estimate
-   `pvalue` : P value, the probability under the assumption of no
    effect or no difference, of obtaining a result equal to or more
    extreme than what was actually observed
-   `padjusted` : A p-value adjustment
-   `logCPM` : log counts per million, measure of expression level <br>


Mobile mRNA Analysis
------------------------------------------------------------------------
For the identification of mobile mRNA, the same `mobileRNA` workflow can be
used. For the alignment step, as well as using a merged genome assembly file
(FASTA) the user will also need to generate a merged genome gene annotation file
(GFF) with the same chromosome labels. As well as a dataframe contain 
sample data where rows represent samples, and column 1 stores the sample names, 
column 2 stores the fastq file name for mate 1 and, for pair-end alignment, 
column 3 stores the fastq file name for mate 2. While for all other downstream 
functions simply change some of the parameters to set the option to mRNA input 
data type. 


Advanced Analysis
------------------------------------------------------------------------
The quick start analysis allows for the retrieval of a data frame the
sRNA total population in the experimental design and also the candidate
mobile sRNAs. Users may want to advance the analysis and plot the data:

-   Exploratory and quality control analysis, such as PCA and distance
    matrices.
-   Summary values including RPM mean and Count mean across specific
    samples.
-   Plotting of the distribution of sRNA classes and the consensus
    dicercall across individual replicates or across the data set.
-   Volcano plot for mRNA data 


Advanced features also include tools to assist functional analysis: 

* Identify genomic features associates with the sRNA clusters (ie. to explore the RNA expression of sRNA-producing genes in parallel analysis)
* Extract consensus RNA sequence for target prediction analysis

Please head to the vignette for the advanced analysis options.
<br> <br>
------------------------------------------------------------------------

*Last updated:* 18-10-2023
