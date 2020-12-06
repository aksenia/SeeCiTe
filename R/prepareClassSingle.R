#' Prepare data summaries for final classifications, a helper function.
#' Performs main steps to check support of CNV in a single sample data
#' @param analyzed_df A dataframe processed with runAnalyseSignal
#'
#' @return A data.table formatted for final classification with classifySingles
#'
#' @examples
prepareClassSingle <- function(analyzed_df){
  # bafType and lrrType Tests are offspring-centered!
  analyzed_df %>%
    dplyr::rename(lrrm = LRR_mean,
                  lrrsd=LRR_SD,
                  flank_lrrsd = flank_LRR_SD_stab,
                  flank_lrrm = flank_LRR_mean_stab) %>%
    dplyr::distinct() %>%
    dplyr::group_by(sample, coordcnv) %>%
    dplyr::summarise(lrrsd_O =  lrrsd[relation=="O"],
                     flank_lrrsd_O =  flank_lrrsd[relation=="O"],
                     lrrm_O =  lrrm[relation=="O"],
                     flank_lrrm_O =  flank_lrrm[relation=="O"],
                     lrrmedtest_O =  lrrmedtest[relation=="O"],
                     BAFprobes_O = BAFprobes[relation=="O"],
                     flanktest = paste(sort(unique(flanktest[cohen_cutoff==TRUE])), collapse=":"),
                     bafTypeTest = unique(bafTypeTest),
                     lrrTypeTest = unique(lrrTypeTest),
                     type=unique(type),
                     copynumber=unique(copynumber),
                     numsnp=unique(numsnp)) %>%
    dplyr::mutate(flanktest = ifelse(flanktest =="", "CNVresemblesflank", flanktest),
                  lrrm_dev_O = ifelse(grepl("Oflank", flanktest),
                                      ifelse(lrrm_O>0, "dir_dup", "dir_del"),
                                      ifelse(grepl("agrees", bafTypeTest) & abs(flank_lrrm_O) > 0.5, "wider_event", NA)))
}
