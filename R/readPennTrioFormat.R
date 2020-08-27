#' Load PennCNV trio adapted file into R. Expects the coding "sample":"relation":"status" in the fifth column.
#'
#' @param filename A path to the output file of PennCNV trio, produced in the infer inheritance steps.
#'
#' @return A table from the input  file
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
readPennTrioFormat <- function(filename) {
  read.table(filename,
                    stringsAsFactors = F,
                    col.names = c("coordcnv", "numsnp", "length", "HMMstate",
                                  "sampleid", "startsnp", "endsnp")) %>%
    tidyr::separate(HMMstate, c("state", "copynumber"), sep = ",") %>%
    tidyr::separate(sampleid, c("sample", "relation", "status"), sep = ":") %>%
    dplyr::mutate(sample = gsub(":+$", "", sample),
                  numsnp = as.integer(gsub(".*=", "", numsnp)),
                  length = as.integer(gsub(".*=|,", "", length)),
                  copynumber = as.integer(gsub(".*=", "", copynumber)),
                  startsnp = gsub(".*=", "", startsnp),
                  endsnp = gsub(".*=", "", endsnp))
}
