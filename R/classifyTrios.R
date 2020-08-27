#' Main function that uses prepared summary statistics to classify CNVs into categories
#' @param analyzed_df Output of runAnalyzeSignal
#' @param lrr_flank_pass A cutoff on the LRR_SD in flanking region that defines placement of putative de novo CNVs into Unlikely category for more detailed examination
#'
#' @return A data table with most useful statistics per individual in a trio and classification results for an offspring CNV
#'
#' @examples
classifyTrios <- function(analyzed_df, lrr_flank_pass=0.2){
  prepare_dt <- prepareClass(analyzed_df)
  class_dt <- prepare_dt  %>%
    dplyr::mutate(Chr = gsub(":.*|chr", "", coordcnv)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(offspring_consistType = ifelse((bafTypeTest!="BAFagreesType") + (lrrTypeTest!="LRRagreesType") > 0, "OconflictCNVType", "OconsistentCNVType"),
                  flanktestShort = ifelse(!grepl("resembles", flanktest),
                                          paste(sort(sapply(strsplit(flanktest, split = ":")[[1]],
                                                            function(s) substr(s,1,1))), collapse = ""),
                                          paste0(c("F","M","O")[!grepl("cnv_sim", c(lrrmedtest_F, lrrmedtest_M, lrrmedtest_O))], collapse="")),
                  flanktestShort = ifelse(flanktestShort == "", "none", flanktestShort),
                  corr_flanktestShort = ifelse(grepl("M", flanktestShort),
                                               ifelse(grepl("agreeNormal", BAFprobes_M) || (lrrmedtest_M == "cnvmd_at_zero" && abs(lrrm_M)<0.1), gsub("M", "", flanktestShort), flanktestShort),
                                               flanktestShort),
                  corr_flanktestShort = ifelse(grepl("F", corr_flanktestShort),
                                               ifelse(grepl("agreeNormal", BAFprobes_F) || (lrrmedtest_F == "cnvmd_at_zero" && abs(lrrm_F)<0.1), gsub("F", "", corr_flanktestShort), corr_flanktestShort),
                                               corr_flanktestShort)) %>%
    dplyr::mutate(inheritanceTest = ifelse(lrrtest == corr_flanktestShort ||
                                             (lrrtest == "FM" & grepl(gsub("O", "", flanktestShort), gsub("O", "", baftest)) & flanktestShort %in% c("FO", "MO") & n_consistentP==1) ||
                                             (lrrtest == flanktestShort & grepl(gsub("O", "", lrrtest), gsub("O", "", baftest)) & offspring_consistType == "OconsistentCNVType") & n_consistentP==1,
                                           ifelse(grepl("mosaic", mosaic_M)||grepl("mosaic", mosaic_F),"mosaic", "inherited"),
                                           ifelse((lrrtest == "FM" & corr_flanktestShort == "O" & n_consistentP==2)||
                                                    (flanktestShort=="O" & n_consistentP==2 & (!grepl(gsub("O", "", lrrtest), gsub("O", "", baftest)) ||
                                                                                                 (baftest == "FM" & grepl("agreeNormal", BAFprobes_M) ))), "denovo",
                                                  ifelse(grepl("mosaic", mosaic_M)||grepl("mosaic", mosaic_F),"mosaic", "unclear"))),
                  inheritanceTest = ifelse(type != "dup" || inheritanceTest != "unclear", inheritanceTest,
                                           ifelse(((baftest == "FO" & grepl("suggDup", BAFprobes_F) &
                                                      (grepl("cnv_gt_flank", lrrmedtest_F) + grepl("cnv_gt_flank", lrrmedtest_O) >=1)) ||
                                                     ((baftest == "MO" & grepl("suggDup", BAFprobes_M) &
                                                         (grepl("cnv_gt_flank", lrrmedtest_M) + grepl("cnv_gt_flank", lrrmedtest_O) >=1)))), "inherited", inheritanceTest)),
                  offspring_flank_lrrsd_pass = ifelse(flank_lrrsd_O <= lrr_flank_pass || status == "inherited" || (status == "ambiguous" & inheritanceTest == "inherited"), "FLANKOK", "FLANKNOISY")) %>%
    tidyr::unite(param, status, inheritanceTest, sep = "_", remove = F) %>%
    dplyr::mutate(param = ifelse(offspring_consistType == "OconflictCNVType" & inheritanceTest!="inherited", "OconflictCNVType", param),
                  param = ifelse(offspring_flank_lrrsd_pass == "FLANKNOISY", "FLANKNOISY", param),
                  seecite=ifelse(grepl("inherited_inherited|denovo_denovo", param), "probable",
                                     ifelse(grepl("FLANKNOISY", param), "unlikely", "borderline")))


}
