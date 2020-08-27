#'  Remove outliers in the CNV locus LRR values, a helper function.
#' @param ds A vector of LRR values
#'
#' @return A vector of regularized LRR values
#'
#' @examples
stabilizeDistr <- function(ds){
  ds <- ds[! is.na(ds)]
  y <- ifelse(length(boxplot.stats(ds)$out)==0, outliers::outlier(ds), boxplot.stats(ds)$out)
  unique(ds[! ds %in% y])
}
