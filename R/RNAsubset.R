#' Subset data to select a specific small RNA population
#'
#' @description Creates a data-frame containing the desired sRNA class(s).
#' See [RNAconsensus()] for information on defining the sRNA type
#' for each cluster. The function allows the choice to filtered the data by
#' statistical significance. Set \code{sig=TRUE} to filtered by significance
#' (p-adjusted).
#' @param data A numerical data-frame containing the sample data and consensus
#' sRNA class/type (see [RNAconsensus()].
#' @param type A number to represent the type of small RNA population to subset
#' for.
#' @param ... Related to number in the `type` argument
#' This can be a value from 20-24. To select, 24-nt sRNA, state 24.
#' Multiple values can be inputted, for instance both 21 and 22 can be
#' stated to select both.
#' @param sig Parameter to filter and select significant sRNA. If
#'  \code{sig=TRUE}, data will be filtered based on p-adjusted < 0.05
#'  significance threshold.
#' @examples
#' data("sRNA_data")
#'
#'##  define consensus sRNA classes.
#'samples <- c("TomEgg_1", "TomEgg_2", "TomEgg_3")
#'
#' # Define consensus
#'sRNA_data_summary <- RNAconsensus(data = sRNA_data,
#'                                      conditions = samples)
#'
#' # Subset data for  24-nt sRNAs
#' sRNA_24 <- RNAsubset(sRNA_data_summary, type = 24)
#'
#'
#'# Subset data for 24 21/22-nt sRNAs
#'sRNA_2122 <- RNAsubset(sRNA_data_summary, type = c(21, 22))
#'
#' # You can subset by any combination of classes. For example, a dataset
#' # of 23-nt & 24-nt sRNAs or just 20-nt sRNAs.
#'
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom dplyr "filter"

RNAsubset <- function(data, type,  sig=FALSE, ...){
  type <- as.character(type)

  if(length(type) == 1){
    if (type == "20"){
      x <- data %>%
        dplyr::filter(sRNA_Consensus == "nt_20")
    } else
      if (type == "21"){
        x <- data %>%
          dplyr::filter(sRNA_Consensus == "nt_21")
      }else
        if (type == "22"){
          x <- data %>%
            dplyr::filter(sRNA_Consensus == "nt_22")
        }else
          if (type == "23"){
            x <- data %>%
              dplyr::filter(sRNA_Consensus == "nt_23")
          }else
            if (type == "24"){
              x <- data %>%
                dplyr::filter(sRNA_Consensus == "nt_24")
            }
  }

  if (length(type) == 2){
    if (.match_vec("20", type) & .match_vec("21", type)) {
      x <- data %>%
        dplyr::filter(sRNA_Consensus == "nt_20" |
                        sRNA_Consensus == "nt_21")
    }else
      if (.match_vec("20", type) & .match_vec("22", type)) {
        x <- data %>%
          dplyr::filter(sRNA_Consensus == "nt_20" |
                          sRNA_Consensus == "nt_22")
      }else
        if (.match_vec("20", type) & .match_vec("23", type)) {
          x <- data %>%
            dplyr::filter(sRNA_Consensus == "nt_20" |
                            sRNA_Consensus == "nt_23")
        }else
          if (.match_vec("20", type) & .match_vec("24", type)) {
            x <- data %>%
              dplyr::filter(sRNA_Consensus == "nt_20" |
                              sRNA_Consensus == "nt_24")
          }else
            if (.match_vec("21", type) & .match_vec("22", type)) {
              x <- data %>%
                dplyr::filter(sRNA_Consensus == "nt_21" |
                                sRNA_Consensus == "nt_22")
            }else
              if (.match_vec("21", type) & .match_vec("23", type)) {
                x <- data %>%
                  dplyr::filter(sRNA_Consensus == "nt_21" |
                                  sRNA_Consensus == "nt_23")
              }else
                if (.match_vec("21", type) & .match_vec("24", type)) {
                  x <- data %>%
                    dplyr::filter(sRNA_Consensus == "nt_21" |
                                    sRNA_Consensus == "nt_24")
                }else
                  if (.match_vec("22", type) & .match_vec("23", type)) {
                    x <- data %>%
                      dplyr::filter(sRNA_Consensus == "nt_22" |
                                      sRNA_Consensus == "nt_23")
                  }else
                    if (.match_vec("22", type) & .match_vec("24", type)) {
                      x <- data %>%
                        dplyr::filter(sRNA_Consensus == "nt_22" |
                                        sRNA_Consensus == "nt_24")
                    }
  }
  if (length(type) == 3){
    if (.match_vec("20", type) & .match_vec("21", type) &
        .match_vec("22", type)) {
      x <- data %>%
        dplyr::filter(sRNA_Consensus == "nt_20" |
                        sRNA_Consensus == "nt_21"|
                        sRNA_Consensus == "nt_22")
    }else
      if (.match_vec("20", type) & .match_vec("21", type) &
          .match_vec("23", type)) {
        x <- data %>%
          dplyr::filter(sRNA_Consensus == "nt_20" |
                          sRNA_Consensus == "nt_21"|
                          sRNA_Consensus == "nt_23")
      }else
        if (.match_vec("20", type) & .match_vec("21", type) &
            .match_vec("24", type)) {
          x <- data %>%
            dplyr::filter(sRNA_Consensus == "nt_20" |
                            sRNA_Consensus == "nt_21"|
                            sRNA_Consensus == "nt_24")
        }else
          if (.match_vec("21", type) & .match_vec("22", type) &
              .match_vec("23", type)) {
            x <- data %>%
              dplyr::filter(sRNA_Consensus == "nt_21" |
                              sRNA_Consensus == "nt_22"|
                              sRNA_Consensus == "nt_23")
          }else
            if (.match_vec("21", type) & .match_vec("22", type) &
                .match_vec("24", type)) {
              x <- data %>%
                dplyr::filter(sRNA_Consensus == "nt_21" |
                                sRNA_Consensus == "nt_22"|
                                sRNA_Consensus == "nt_24")
            }else
              if (.match_vec("22", type) & .match_vec("23", type) &
                  .match_vec("24", type)) {
                x <- data %>%
                  dplyr::filter(sRNA_Consensus == "nt_22" |
                                  sRNA_Consensus == "nt_23"|
                                  sRNA_Consensus == "nt_24")
              }
  }
  if(length(type) == 4){
    if (.match_vec("20", type) & .match_vec("21", type) &
        .match_vec("22", type)& .match_vec("23", type)) {
      x <- data %>%
        dplyr::filter(sRNA_Consensus == "nt_20" |
                        sRNA_Consensus == "nt_21" |
                        sRNA_Consensus == "nt_22"|
                        sRNA_Consensus == "nt_23")
    } else
      if (.match_vec("20", type) & .match_vec("21", type) &
          .match_vec("22", type)& .match_vec("24", type)) {
        x <- data %>%
          dplyr::filter(sRNA_Consensus == "nt_20" |
                          sRNA_Consensus == "nt_21" |
                          sRNA_Consensus == "nt_22"|
                          sRNA_Consensus == "nt_24")
      } else
        if (.match_vec("20", type) & .match_vec("21", type) &
            .match_vec("23", type)& .match_vec("24", type)) {
          x <- data %>%
            dplyr::filter(sRNA_Consensus == "nt_20" |
                            sRNA_Consensus == "nt_21" |
                            sRNA_Consensus == "nt_23"|
                            sRNA_Consensus == "nt_24")
        } else
          if (.match_vec("20", type) & .match_vec("22", type) &
              .match_vec("23", type)& .match_vec("24", type)) {
            x <- data %>%
              dplyr::filter(sRNA_Consensus == "nt_20" |
                              sRNA_Consensus == "nt_22" |
                              sRNA_Consensus == "nt_23"|
                              sRNA_Consensus == "nt_24")
          } else
            if (.match_vec("21", type) & .match_vec("22", type) &
                .match_vec("23", type)& .match_vec("24", type)) {
              x <- data %>%
                dplyr::filter(sRNA_Consensus == "nt_21" |
                                sRNA_Consensus == "nt_22" |
                                sRNA_Consensus == "nt_23"|
                                sRNA_Consensus == "nt_24")
            }
  }
  if (length(type) == 5){
    if (.match_vec("20", type) & .match_vec("21", type) &
        .match_vec("22", type) & .match_vec("23", type) &
        .match_vec("24", type)) {
      x <- data %>%
        dplyr::filter(sRNA_Consensus == "nt_20" |
                        sRNA_Consensus == "nt_21" |
                        sRNA_Consensus == "nt_22"|
                        sRNA_Consensus == "nt_23"|
                        sRNA_Consensus == "nt_24")
    }
  }

  if(sig){
    x<- x %>%
      dplyr::filter(padj < 0.05)
  }
  return(x)
}

