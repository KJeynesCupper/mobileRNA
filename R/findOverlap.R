#' Identify genomic annotations which overlap with the coordinates of the
#' dicer-derived clusters.
#'
#' @description A function to call an OS command to use `BEDtools`
#' to identify the intersect between the merged reference annotation and
#' the merged dicer-derived cluster information in a `.gff3 `format.
#' @param annotation Local directory path to merged reference annotation file
#' (`.gff`).
#' @param out Local directory path to the location to store the file which will
#' be named "genomic_overlap.bed"
#'
#' @param clusters `.gff3` annotation file containing information of the
#' siRNA dicer-derived clusters. This file will have been produced as a part of
#' the pre-processing steps. Specifically, using the output from the first
#' mapping step.
#' @return A file called `"genomic_overlap.bed"` containing annotation
#' information on each dicer-derived cluster.
#' @examples
#' \dontrun{
#' Anno <- "./annotation/merge/anno_merged.gff"
#' save <-   "./cluster_genomic_locations/"
#' ClusterInfo <- "./clusters/ClustersInfo.gff3"
#' findOverlap(annotation = Anno, clusters = ClusterInfo, out = save)
#' }
#' @export
findOverlap <- function(annotation, clusters, out) {
  dir.create(file.path(getwd(),"/cluster_genomic_locations"),
             showWarnings = FALSE)
  system(paste0("bedtools intersect -wa -wb -a", clusters, "-b ",
                annotation," >" , out, "genomic_overlap.bed"  ))
}



