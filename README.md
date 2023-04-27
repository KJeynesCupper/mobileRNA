RNAlocate 
======================================================================

Overview
--------

RNAlocate is an `R` package that provides a pipeline for the rapid identification of mobile RNA molecules in 
plant graft systems. The tool provides an improved pipeline for pre-processing and analysis of sequencing data. 

Author
--------
Katie Jeynes-Cupper, University of Birmingham, kej031@student.bham.ac.uk


Table of Contents
--------
- [Installation](#installation)
- [Loading test data](#Loading-test-data)
- [Getting help](#Getting-help)
- [Pre-mapping](#Pre-mapping)
- [Mapping](#Mapping)
- [Post-mapping analysis](#Post-mapping-analysis)
- [Output](#Output)


Installation
------------

The latest version of the package can be install directly from this GitHub repo:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("KJeynesCupper/RNAlocate")
```

Loading test data
-----------------

To simulate the usage of the package, there is a grafting data set, from a 
tomato-eggplant grafting experiment. The data has already undertaken the 
pre-mapping and mapping steps (pre-processing) and can be found here
`sRNA_data`. The data set has already been organised by the `RNAimport()` 
function available in the package. 

To load the data set, use the following command:

``` r
data("sRNA_data")
```


Getting help
------------

For additional information on each function, please read through the 
documentation in the `RNAlocate` package by typing the `?` help operator before 
any of the function names in the package or by using the `help()` function.

For a step-by-step quick start analysis, consider reading the vignette provided 
with this package:

``` r
vignette("RNAlocate-Quick-Start")
```

For a more in-depth analysis, consider reading the full vignette provided with 
this package:


``` r
vignette("RNAlocate")
```

Pre-mapping
--------------------------------------------
Raw fastq files should be trimmed to remove adapter sequences and low quality 
reads as per best practice. 

#### Installation of Linux Dependencies
The majority of the pre-processing functions require a system call 
to Linux software. The following dependencies need to be installed before you 
begin:

- `FastQC` (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- `bedtools` (https://bedtools.readthedocs.io/en/latest/)
- `Samtools` (http://www.htslib.org/)
- `sed` (https://www.gnu.org/software/sed/manual/sed.html)
- `ShortStack` (https://github.com/MikeAxtell/ShortStack)


Mapping
--------------------------------------------
Here, we introduce an alternative mapping method for the analysis of plant
heterograft samples. The heterograft system involves two genotypes; here
the two genome references are merged into a single reference to which 
samples are aligned to. 

`RNAlocate` offers a function to merge two reference genomes into one: 

``` r
RNAlocate::mergeFiles(files = "/Users/user1/projectname/workplace/reference/*",
        out = "/Users/user1/projectname/workplace/reference/merge/merged_reference.fa")
```

The analysis pipeline is formulated to analyse mapping and clustering results 
produced by `ShortStack` (https://github.com/MikeAxtell/ShortStack). 

To distinguish between the reference genomes in a merged file, it is important 
to make sure the chromosome names between the genomes are different and 
distinguishable. 

Here, we recommend a double-mapping process using `ShortStack`, the steps are 
as follow: 

#### Step 1 - Identify loci of dicer-derived sRNA cluster in each sample 

``` bash
ShortStack \
--readfile <control_1.fastq> \
--genomefile <merged_reference.fa> \
--bowtie_cores 6 \
--mmap n \
--mismatches 0 \
--nohp \
--outdir <./output/directory>

```
#### Step 2 - Uniquely map samples to all loci of identified sRNA clusters

From the output of Step 1, the `identifyClusters()` function can collate all 
the loci information into a single `.gff3` file.  

``` r
sample_names <- c("<treatment_1>", "<treatment_2>", "<control_1>","<control_2>")

folder <- <./output/directory/from/step/1/>
save_folder <- <./output/directory/ClustersInfo.gff3>


RNAlocate::identifyClusters(files = folder, 
             out = save_folder,
             samples = sample_names)
```

Each sample is mapped to the merged reference genome with the annotation file 
containing the cluster loci to analyse. 

``` bash
ShortStack \
--readfile <control_1.fastq> \
--genomefile <merged_reference.fa> \
--locifile <./output/directory/ClustersInfo.gff3> \
--bowtie_cores 6 \
--mmap n \
--mismatches 0 \
--nohp \
--mincov 5 \
--outdir <./output/directory/step2/>

```

 
Post-mapping analysis 
--------------------------------------------
The aligned data can now be analysed in R, and potential mobile sRNA can be 
identified. 


#### Step 1: Import data 
State the location of the mapping results, the sample names and sRNA cluster loci
annotation file. 

``` r
# Directory containing results
results_dir <-  "<./output/directory/step2/>"

# Sample names and total number of reads, in the same order. 
sample_names <- c("<treatment_1>", "<treatment_2>", "<control_1>","<control_2>")

# sRNA cluster loci annotation file 
clusters_anno <- rtracklayer::import.gff(<./output/directory/ClustersInfo.gff3>)


sRNA_data <- RNAimport(results = results_dir, 
                           samples = sample_names, 
                           clusters = clusters_anno)
                           
```


#### Step 2: Calculate the consensus of each sRNA cluster  
Each sample in the analysis has determined the class of each 
dicer-derived sRNA cluster based on the most abundant small RNA size in the 
sample. It is expected that replicates within the same condition will define the 
same cluster as the same class. The `RNAconsensus()` function is used to define 
the class of a sRNA cluster based on the consensus across specific replicates. 

To identify potentially mobile sRNA moving from genotype A to B in a heterograft,
in comparison to a B-genotype self-graft, it is recommended to base 
the consensus call on the heterograft samples. This will ensure that the 
sRNA class is more accurately defined by the genotype it originates from. 

Below, replicates "treatment_1" and "treatment_2" represent two heterograft 
samples, while "control_1" and "control_2" represent self-graft samples.

``` r

samples <- c("<treatment_1>", "<treatment_2>")

sRNA_data_summary <- RNAconsensus(data = sRNA_data, 
                                 conditions = samples, 
                                 tidy=TRUE)

```


#### Step 3: Identify potential mobile sRNA 
Finally, the `RNAmobile()` function filters the dataset to retain the sRNA 
clusters mapped to the A-genome. It is important, that each genome is 
distinguishable by the chromosome names, hence, in this example the A-genome
contains "A" before the chromosome number. 

``` r
# define control samples
controls <- c("<control_1>", "<control_2>")

mobile <- RNAmobile(data = sRNA_data_summary, 
                    controls = controls,
                    id = "A", 
                    task = "keep")

```


Output 
--------------------------------------------
A dataframe where rows represent potential mobile sRNA clusters.
The columns include information on the cluster, individual sample replicates, 
and more. 

### Information on the cluster:
- `chr`: Name of the chromosome or scaffold
- `start` : Start position of the cluster
- `end` : End position of the cluster
- `width` : Length of the locus (base-pairs)

### Information on each sample replicate:
- `cluster`: Cluster Name
- `Dicer-call` : The size of most abundant small RNA size
- `Count` : Number of uniquely aligned (*e.g.* not multi-mapping) reads that 
overlap this locus.
- `RPM` : Reads per Million


### Other information
- `sRNA_Consensus` : Consensus sRNA class calculated by `RNAconsensus()`

If an annotation file was imported and overlapped using `RNAimport()`, parameter
- `source` : name of the program that generated this feature, or the data source 
(database or project name)
- `feature` : feature type name, e.g. Gene, Variation, Similarity
- `score` : A floating point value
- `strand` : defined as + (forward) or - (reverse)
- `frame` : One of '0', '1' or '2'. '0' indicates that the first base of the 
feature is the first base of a codon, '1' that the second base is the first 
base of a codon, and so on
- `attribute` : A semicolon-separated list of tag-value pairs, providing 
additional information about each feature.

If differential expression analysis was undertaken with`RNAanalysis()`: 
- `log2FoldChange` : Log2FoldChange–The effect size estimate
- `pvalue` :  P value, the probability under the assumption of no effect or 
no difference, of obtaining a result equal to or more extreme than what was 
actually observed
- `padjusted` : A p-value adjustment
- `logCPM` : log counts per million, measure of expression level

If the mean RPM and Count was calculated `RNAmean()`: 
- `mean_RPM` : mean RPM, based on parameters 
- `mean_Count` : mean counts, based on parameters 
------------------------------------------------------------------------

*Last updated:* 26-04-2023
