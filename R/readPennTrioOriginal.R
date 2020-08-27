#' Load PennCNV trio file in original format as supplied by PennCNV into R
#'
#' @param filename A path to the output file of PennCNV trio
#' @param triostate Is a column with trio state present? Default TRUE.
#'
#' @return A table with the file content.
#' @export
#'
#' @importFrom magrittr "%>%"
#' @examples
readPennTrioOriginal <- function(filename, triostate=TRUE) {
  if(triostate==TRUE){
    header <- c("coordcnv", "numsnp", "length", "HMMstate",
                                "sample", "startsnp", "endsnp", "relation", "triostate")
  } else {
    header <- c("coordcnv", "numsnp", "length", "HMMstate", "sample", "startsnp", "endsnp")
}
  result <- read.table(filename,
                    stringsAsFactors = F,
                    col.names = header) %>%
    tidyr::separate(HMMstate, c("state", "copynumber"), sep = ",") %>%
    dplyr::mutate(sample=basename(sample),
                  numsnp = as.integer(gsub(".*=", "", numsnp)),
                  length = as.integer(gsub(".*=|,", "", length)),
                  copynumber = as.integer(gsub(".*=", "", copynumber)),
                  startsnp = gsub(".*=", "", startsnp),
                  endsnp = gsub(".*=", "", endsnp)) %>%
    tidyr::separate(coordcnv, c("chr", "coord"), sep = ":", remove = F) %>%
    tidyr::separate(coord, c("start", "end"), sep = "-", remove = T) %>%
    dplyr::mutate(start = as.integer(start),
           end = as.integer(end),
           length = as.integer(length))
  if (triostate==TRUE) {
    result <- result %>%
      dplyr::mutate(triostate = gsub(".*=", "", triostate)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(triostate = paste(sort(unique(strsplit(triostate, split="-")[[1]])), collapse="-"))
  }
  result
}
