#' Prepare data summaries for final classifications, a helper function.
#' Performs main steps to check support of CNV in offspring and parent data
#' @param analyzed_df A dataframe processed with runAnalyseSignal
#'
#' @return A data.table formatted for final classification with classifyTrios
#'
#' @examples
prepareClass <- function(analyzed_df){
  # bafType and lrrType Tests are offspring-centered!
  analyzed_df %>%
    dplyr::rename(lrrm = LRR_mean,
                  lrrsd=LRR_SD,
                  flank_lrrsd = flank_LRR_SD_stab,
                  flank_lrrm = flank_LRR_mean_stab) %>%
    dplyr::distinct() %>%
    dplyr::mutate(consistentParent =
                    ifelse(cohen_cutoff==FALSE,
                           ifelse(sum(BAFprobes %in% c("BAFagreeNormal", "possiblyLOH"))>0, T, F),
                           ifelse(BAFprobes == "BAFagreeNormal" || lrrmedtest == "cnvmd_at_zero", T, F))) %>%
    dplyr::group_by(sample, coordcnv) %>%
    dplyr::summarise(lrrsd_O =  lrrsd[relation=="O"],
                     lrrsd_F =  lrrsd[relation=="F"],
                     lrrsd_M =  lrrsd[relation=="M"],
                     flank_lrrsd_O =  flank_lrrsd[relation=="O"],
                     flank_lrrsd_F =  flank_lrrsd[relation=="F"],
                     flank_lrrsd_M =  flank_lrrsd[relation=="M"],
                     lrrm_O =  lrrm[relation=="O"],
                     lrrm_F =  lrrm[relation=="F"],
                     lrrm_M =  lrrm[relation=="M"],
                     flank_lrrm_O =  flank_lrrm[relation=="O"],
                     flank_lrrm_F =  flank_lrrm[relation=="F"],
                     flank_lrrm_M =  flank_lrrm[relation=="M"],
                     lrrmedtest_O =  lrrmedtest[relation=="O"],
                     lrrmedtest_F =  lrrmedtest[relation=="F"],
                     lrrmedtest_M =  lrrmedtest[relation=="M"],
                     BAFprobes_O = BAFprobes[relation=="O"],
                     BAFprobes_M = BAFprobes[relation=="M"],
                     BAFprobes_F = BAFprobes[relation=="F"],
                     n_consistentP = sum(consistentParent[relation!="O"]),
                     flanktest = paste(sort(unique(flanktest[cohen_cutoff==TRUE])), collapse=":"),
                     baftest = unique(baftest),
                     lrrtest = unique(lrrtest),
                     bafTypeTest = unique(bafTypeTest),
                     lrrTypeTest = unique(lrrTypeTest),
                     type=unique(type),
                     copynumber=unique(copynumber),
                     parental_orig=unique(parental_orig),
                     numsnp=unique(numsnp),
                     length=unique(length),
                     status=unique(status),
                     countBafTypeDisagrP = sum(!grepl("agreeNormal|LOH", unique(BAFprobes[relation!="O"]))),
                     bafTypeDisagrP = paste(unique(BAFprobes[!grepl("agreeNormal|LOH", BAFprobes[relation!="O"])]), collapse=", ")) %>%
    dplyr::mutate(bafTypeDisagr = ifelse(countBafTypeDisagrP == 0, "BAFTypeAgreeP", bafTypeDisagrP),
                  flanktest = ifelse(flanktest =="", "OFM:CNVresemblesflank", flanktest),
                  lrrm_dev_O = ifelse(grepl("Oflank", flanktest),
                                      ifelse(lrrm_O>0, "dir_dup", "dir_del"),
                                      ifelse(grepl("agrees", bafTypeTest) & abs(flank_lrrm_O) > 0.5, "wider_event", NA)),
                  lrrm_dev_M = ifelse(grepl("Mflank", flanktest),
                                      ifelse(lrrm_M>0, "dir_dup", "dir_del"),
                                      ifelse(abs(flank_lrrm_M) > 0.2,
                                             ifelse(lrrm_M>0, "dir_dup", "dir_del"),
                                             "no_dev")),
                  lrrm_dev_F = ifelse(grepl("Fflank", flanktest),
                                      ifelse(lrrm_F>0, "dir_dup", "dir_del"),
                                      ifelse(abs(flank_lrrm_F) > 0.2,
                                             ifelse(lrrm_F>0, "dir_dup", "dir_del"),
                                             "no_dev")),
                  mosaic_M = ifelse(lrrm_dev_M=="dir_del" & grepl("BAFsuggDup", BAFprobes_M), "mosaic", NA),
                  mosaic_F = ifelse(lrrm_dev_F=="dir_del" & grepl("BAFsuggDup", BAFprobes_F), "mosaic", NA)) %>%
    dplyr::select(-bafTypeDisagrP)
}
