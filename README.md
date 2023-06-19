mobileRNA <a href="kjeynescupper.github.io/mobileRNA/"><img src="man/figures/logo.png" align="right" height="138" /></a>
======================================================================
<br>

Overview
======================================================================
<br>
mobileRNA is an `R` package that provides a pipeline for the rapid 
identification of endogenous mobile small RNA (sRNA) molecules in plant graft 
systems. The tool provides a pipeline for pre-processing and analysis of 
sRNA sequencing data, and soon mRNA sequencing data. 

It has been established that many different substances and molecules 
including RNAs can travel across the graft junction. Plant heterograft systems are 
comprised of two genotypes joined at the graft junction; hence, molecules produces
and encoded by each genotype can move across the graft junction and be exchanged. 
These molecules could have implications to the regulation of gene expression and 
traitacquisition. 

Current methods utilise a step-wise mapping of samples to each genome within the
graft system. While, here we introduce a new mapping method where we align 
each sample replicates to a merge genome reference comprised of both genome 
assemblies relating to the genotypes in the heterograft system. 

**Look-out for Version 2 which accommodates mRNA movement**

Author
--------
Katie Jeynes-Cupper, University of Birmingham, kej031@student.bham.ac.uk


Table of Contents
--------
- [Installation](#installation)
- [Loading Test Data](#Loading-test-data)
- [Getting Help](#Getting-Help)
- [Summary](#Summary)
- [Merging Genome Assemblies](#Merging-Genome-Assemblies)
- [Auto-Detection of sRNA Cluster](#Auto-Detection-of-sRNA-Cluster)
- [Mapping](#Mapping)
- [Analysis](#Analysis)
- [Output](#Output)
- [Functional Analysis](#Functional-Analysis)
- [Optional Extras](#Optional-Extras)

<br>
Installation
------------

The latest version of the package can be install directly from this GitHub repo:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("KJeynesCupper/mobileRNA", branch = "main")
```
<br>
Loading Test Data
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

<br>

Getting Help
------------

For additional information on each function, please read through the 
documentation in the `mobileRNA` package by typing the `?` help operator before 
any of the function names in the package or by using the `help()` function.

For an in-depth step-by-step analysis, consider reading the vignette provided 
with this package:


``` r
vignette("mobileRNA")
```
<br>

Summary
------------
<p>
    <img src="./man/figures/program_flow.png" width="400" height="350" align="right" />
</p>

The workflow is shown in the figure to the right. 
It begin in R-Studio to merge the two genome assemblies into one, then the 
pre-processing moves into Linux to align each replicate to the merged reference 
and then back into R-Studio to undertake the analysis to identify potentially 
mobile RNA species.  

Going forward,we assume standard quality control steps on raw samples has been
completed (i.e. trimming of adapters and low quality reads)
  

<br>
<br>
<br>
<br>

Merging Genome Assemblies
--------------------------------------------
Our method is built on aligning samples to a merged reference genome, comprised 
of both genome assemblies relating to the genotypes in the heterograft system. 
Merge the two genomes assemblies into a single file using the function in R 
below. 

``` r
merged_reference <- RNAmergeGenomes(genomeA = "./workplace/reference/ref1.fa",
                                    genomeB = "./workplace/reference/ref2.fa",
                        out_dir = "./workplace/reference/merge/merged_ref.fa")

```


<br>

Auto-Detection of sRNA Cluster  
--------------------------------------------
<br>

Here we identify and build a list sRNA clusters within each sample to assist
the mapping step later on to ensure consistency across the analysis.

We recommend installing the `ShortStack` 
(https://github.com/MikeAxtell/ShortStack) program to detect sRNA clusters and 
align sRNAseq samples. 



#### Step 1 - Cluster analysis with ShortStack

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
<br>

#### Step 2 - Build sRNA cluster list 

Now, we collate all the sRNA loci information from each sample into a text file. 

``` r
# location of step 1 output
folder <- <./output/directory/from/step/1/>

# name and location to save output file to (must be .txt)
save_folder <- <./output/directory/ClustersInfo.txt>

# names of samples (ie. folder names)
sample_names <- c("<treatment_1>", "<treatment_2>", "<control_1>","<control_2>")


loci_info <- RNAloci(files = folder, 
             out = save_folder,
             samples = sample_names)
```
<br>

Mapping   
--------------------------------------------
Each sample is mapped to the merged reference genome with the list of sRNA 
clusters. 

<br>

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

<br>

Analysis 
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


sRNA_data <- RNAimport(input = "sRNA", 
                       directory = results_dir,
                       samples = sample_names)
                           
```

<br>

#### Step 2: Calculate the consensus of each sRNA cluster  
For a given sRNA cluster, each replicate has determined the class (20-24nt) 
based on the most abundant small RNA size. Replicates within the same condition 
are expected to class a given sRNA similarly. 

The `RNAconsensus()` function is used to define the class of a sRNA cluster 
based on the consensus across specific replicates. To identify forigen mobile 
sRNAs, it is recommended to base the consensus call on the heterograft samples. 
This will ensure that the sRNA class is more accurately defined by the genotype 
it originates from. 

Below, replicates "treatment_1" and "treatment_2" represent two heterograft 
samples, while "control_1" and "control_2" represent self-graft samples.

``` r

samples <- c("<treatment_1>", "<treatment_2>")

sRNA_data_summary <- RNAconsensus(data = sRNA_data, 
                                 conditions = samples, 
                                 tidy=TRUE)

```
<br>

#### Step 3: Identify potential mobile sRNA 
Finally, the `RNAmobile()` function filters the dataset to retain the sRNA 
clusters mapped to the A-genome. It is important, that each genome is 
distinguishable by the chromosome names, hence, in this example the A-genome
contains "A" before the chromosome number. 

``` r
# define control samples
controls <- c("<control_1>", "<control_2>")

mobile_df <- RNAmobile(data = sRNA_data_summary, 
                    controls = controls,
                    id = "A", 
                    task = "keep")

# output dataframe containing potentially mobile sRNAs
output <- mobile_df

# save output as txt file 
write.table(output, "./output.txt")


```
<br>

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

<br>
### Other information
- `sRNA_Consensus` : Consensus sRNA class calculated by `RNAconsensus()`

If an annotation file was imported and overlapped using `RNAattributes()`:

- `source` : name of the program that generated this feature, or the data source (database or project name)
- `feature` : feature type name, e.g. Gene, Variation, Similarity
- `score` : A floating point value
- `strand` : defined as + (forward) or - (reverse)
- `frame` : One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on
- `attribute` : A semicolon-separated list of tag-value pairs, providing additional information about each feature.

If differential expression analysis was undertaken with`RNAanalysis()`: 

- `log2FoldChange` : Log2FoldChange–The effect size estimate
- `pvalue` :  P value, the probability under the assumption of no effect or no difference, of obtaining a result equal to or more extreme than what was actually observed
- `padjusted` : A p-value adjustment
- `logCPM` : log counts per million, measure of expression level

If the mean RPM and Count was calculated `RNAmean()`: 
- `mean_RPM` : mean RPM, based on parameters 
- `mean_Count` : mean counts, based on parameters 

<br>

Functional Analysis 
--------------------------------------------
By now you will have retrieved a list of potentially mobile sRNA molecules and 
you will want to identify whether they play a role in the biological system. 

Small interfering RNAs (siRNAs) have a length of 24-nucleotides, and are known
to play a role in the RNA-directed DNA methylation pathway via targeting 
complementary sequencing in the DNA.

Using the genomic location of the sRNA clusters we can identify overlaps in 
the genome of the origin tissue to predict the implication of the 
mobile sRNA. This function returns the input dataframe with additional fields 
if overlaps were located:

``` r
mobile_df_attributes <- RNAattributes(data = mobile_df,
                            annotation = "./annotation/origin_annotation.gff3")

```
<br>

We can also locate specific genomic features which overlap with the the sRNA 
cluster locus including promoter regions, exon, introns, untranslated regions 
and repeat regions:


``` r
mobile_df_features <- RNAfeatures(data = mobile_df,
                            annotation = "./annotation/origin_annotation.gff3", 
                      repeats = "./annotation/origin_annotation_repeats.gff3")

```
<br>


Lastly, lets extract the sRNA cluster sequences.These can utilised with `BLAST`, 
or plant sRNA target prediction tools such as `TargetFinder` 
(https://github.com/carringtonlab/TargetFinder) or overlapped with `mRNAseq` 
data to identify complementary sequences to further elucidate the potential 
function of sRNAs. 

The `RNAsequences` function identifies whether the most abundant sRNA is 
consistent across the replicates, and if so, it extracts the sRNA nucleotide
sequence and calculates the RNA and DNA complementary sequences, as well as 
stating the length of the sequence. 
``` r

mobile_sequences <- RNAsequences(mobile_df)

# save output as txt file 
write.table(output, "./output.txt")

```

<br> 


Optional Extras
--------------------------------------------
The package also includes functions for: 

* Exploratory and quality control analysis, such as PCA & heatmap plots. 
* Summary values including RPM mean and Count mean across specific samples.
* Plot the distribution of sRNA classes (20-24nt) across individual replicates or across the dataset. 
* Statistical analysis using differential methods from either DESeq2 or edgeR. 

The package workflow can easily be manipulated to enable the identification of local populations of RNA species. 

------------------------------------------------------------------------

*Last updated:* 06-06-2023
