#' Main function that uses prepared summary statistics to classify CNVs into categories
#' @param analyzed_df Output of runAnalyzeSignal
#' @param lrr_flank_pass A cutoff on the LRR_SD in flanking region that defines placement of putative de novo CNVs into Unlikely category for more detailed examination
#'
#' @return A data table with most useful statistics per individual in a trio and classification results for an offspring CNV
#' @export
#' @examples
classifySingles <- function(analyzed_df, lrr_flank_pass=0.3){
  prepare_dt <- prepareClassSingle(analyzed_df)
  class_dt <- prepare_dt  %>%
    dplyr::mutate(Chr = gsub(":.*|chr", "", coordcnv)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(offspring_consistType = ifelse((bafTypeTest!="BAFagreesType") + (lrrTypeTest!="LRRagreesType") > 0, "OconflictCNVType", "OconsistentCNVType"),
                  flanktestShort = gsub("O", "", flanktest)) %>%
    dplyr::mutate(offspring_flank_lrrsd_pass = ifelse(flank_lrrsd_O <= lrr_flank_pass, "FLANKOK", "FLANKNOISY"),
                  param = offspring_consistType,
                  param = ifelse(offspring_flank_lrrsd_pass == "FLANKNOISY", "FLANKNOISY", param),
                  seecite=param,
                  param=flanktestShort)

  class_dt
}
