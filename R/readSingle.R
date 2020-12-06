#' Load and format data for single sample processing
#'
#' @param snp_file A file with extracted snp data, after running 'extract_snp_single.py'
#'
#' @return A data.frame suitable for downstream single sample analysis
#' @export
#'
#'
readSingle <- function(snp_file){
  read.table(snp_file,
             stringsAsFactors = F,
             header = T,
             sep="\t") %>%
    dplyr::mutate(relation="O",
                  Origin=NA,
                  type= ifelse(copynumber<2, "del", "dup")) %>%
    tidyr::gather(-Name, -Chr, -Position, -locus, -copynumber, -type, -numsnp, -Origin,
                  -coordcnv, -sample, -relation, key="parameter", value="value") %>%
    dplyr::group_by(Name, Chr, Position, locus, coordcnv, sample, copynumber, numsnp,
                    relation, Origin, type, parameter) %>%
    dplyr::summarise(value=unique(value)) %>%
    dplyr::ungroup()
}


