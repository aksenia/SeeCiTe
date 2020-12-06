#' SeeCiTe analysis of BAF and LRR in candidate offspring CNV regions in a trio
#' @param input_data The input data table collated in "data" argument with readInput for individual CNV and sample.
#' @param cohen_dist A cutoff on Cohen distance between CNV and flank LRR value distributions. If cohen distance is equal or smaller that this cutoff the distributions are considered not different without further testing. Default is 3.
#' @param sign_lrr_dist A cutoff of p-value at which difference between two individual's CNV LRR value distributions is set as significant, using Kernel density based global two-sample comparison test from ks package. Default is 0.01.
#' @param sign_correct_baf A cutoff p-value at which one-sided wilcoxon test of the shift of CNV LRR distribution from zero is considered significant. Used as an additional evidence to correct status based on BAF values which is succeptible to errors when signal is noisy. Default is 0.01.
#' @param percent_opposite A cutoff on percentage of probes that have LRR values in direction contradicting the suggested CNV type, eg positive for a del and negative for a dup. Default is 20.
#'
#' @return A data table with test results (p-values and classifications) as well as various summaries for CNV and flanking regions per trio and CNV
#'
#' @examples
analyzeSignalConcordance <- function(input_data, cohen_dist=3, sign_lrr_dist=0.01, sign_correct_baf=0.01, percent_opposite=20){
  ## BAF clustering
  baf <- input_data %>%
    dplyr::filter(parameter == "BAF", locus == "CNV") %>%
    dplyr::select(Name, value, relation)
  dbaf <- classifyBAF(baf)
  # LRR testing
  blrr <- input_data %>%
    dplyr::filter(parameter == "LRR")

  bdt <- blrr %>%
    dplyr::filter(parameter == "LRR") %>%
    dplyr::select(Name, Position, value, relation, locus) %>%
    tidyr::spread(relation, value) %>%
    tidyr::unite(Uname, Name, Position, sep = ":") %>%
    tidyr::gather(variable, val, -(Uname:locus)) %>%
    tidyr::unite(temp, locus, variable) %>%
    tidyr::spread(temp, val) %>%
    tidyr::separate(Uname, c("Name", "Position"), sep = ":") %>%
    dplyr::rename(Fa = CNV_F,
                  M = CNV_M,
                  O=CNV_O)

  # run tests to compare distributions
  ddtst <- testLRRDistributions(bdt, cohen_dist=cohen_dist, sign_lrr_dist=sign_lrr_dist)

  flanktest <- ddtst %>%
    dplyr::filter(grepl("flank", lrrtest)) %>%
    dplyr::rename(pvalFlank = ks_pval, flanktest = lrrtest) %>%
    dplyr::select(-dplyr::matches("_pval"), -pvalLRR, -test) %>%
    dplyr::mutate(relation = gsub("flank", "", flanktest))

  meanstest <- ddtst %>%
    dplyr::filter(grepl("flank", lrrtest)) %>%
    dplyr::select(lrrtest, dplyr::matches("lrr_wx"))

  dtst <- ddtst %>%
    dplyr::filter(!grepl("flank", lrrtest))
  # least different pair
  spr <- unlist(dtst[nrow(dtst), "lrrtest"], use.names = FALSE)

  dist_spr <- dtst %>%
    dplyr::ungroup() %>%
    dplyr::filter(hellinger == min(hellinger)) %>%
    dplyr::select(lrrtest) %>%
    unlist(use.names=F)

  dist_spr <- paste0(dist_spr, sep = ":")

  # LRR summary
  lrrS <- blrr %>%
    dplyr::group_by(coordcnv, relation, locus) %>%
    dplyr::summarise(LRR_mean = mean(value, na.rm = T),
              LRR_mean_stab = mean(stabilizeDistr(value), na.rm = T),
              LRR_median = median(value, na.rm = T),
              LRR_median_stab = median(stabilizeDistr(value), na.rm = T),
              LRR_SD = sd(value, na.rm = T),
              LRR_SD_stab = sd(stabilizeDistr(value), na.rm = T),
              N_POS_LRR = sum(value > 0, na.rm = T),
              N_POS_LRR_stab = sum(stabilizeDistr(value) > 0, na.rm = T),
              N_NEG_LRR = sum(value < 0, na.rm = T),
              N_NEG_LRR_stab = sum(stabilizeDistr(value) < 0, na.rm = T)) %>%
    dplyr::ungroup() %>%
    tidyr::gather(variable, val, -coordcnv,-relation,-locus) %>%
    tidyr::unite(temp, locus, variable) %>%
    dplyr::mutate(temp = gsub("CNV_", "", temp)) %>%
    tidyr::spread(temp, val)
  # BAF clustering patterns
  # check how many family members have same number
  # of BAF clusters
  # with normal clustering (LOH clistering check - basically one always must have 0 and 1 cluster)
  clpattn <- dbaf %>%
    dplyr::select(relation, BAFprobes) %>%
    dplyr::distinct() %>%
    dplyr::group_by(BAFprobes) %>%
    dplyr::summarise(rel = paste(sort(relation), collapse = ","), n_mem = dplyr::n())

  # In case of coinsiding clustering patterns
  if (sum(clpattn$n_mem > 1) > 0) {
    # pair with coinciding BAF clusters
    clbafset <- sort(strsplit(clpattn %>%
                                dplyr::filter(n_mem > 1) %>%
                                dplyr::select(rel)  %>%
                                unlist(), split = ",")[[1]])
    clbafset <- paste(clbafset, collapse = "")
  } else if (sum(grepl("possiblyLOH", clpattn$BAFprobes))) {
    # LOH alone is ambiguous with respect to copy number change
    clbafset <- "ambiguous"
  } else {
    clbafset <- "none"
  }
  # pair with more similar LRR
  #  cllrrset <- sort(c(substr(spr, 1, 1), substr(spr, 2, 2)))
  # based on Hellinger distance
  cllrrset <- sort(c(substr(dist_spr, 1, 1), substr(dist_spr, 2, 2)))
  cllrrset <- paste(cllrrset, collapse = "")
  # create BAF
  lrrstrDF <- dtst %>%
    dplyr::mutate(str = "LRRtest",
           str1 = "LRRdist") %>%
    dplyr::mutate(pvalLRRs = format(pvalLRR, scientific = T, digits = 2),
           hellingerStr = format(hellinger, scientific = T, digits = 2)) %>%
    tidyr::unite(lrr_str, lrrtest, pvalLRRs, sep = ":", remove = F) %>%
    tidyr::unite(lrrdist_str, lrrtest, hellingerStr, sep = ":", remove = F)

  lrrStr <- paste(unique(lrrstrDF$test), "test", paste0(lrrstrDF$lrr_str, collapse  = " "))

  lrrDStr <- with(lrrstrDF %>% dplyr::arrange(-hellinger), paste(unique(str1),
                                                          paste0(lrrdist_str, collapse  = " ")))
  lrrMeans <- paste("LRR mean",
                    paste0((lrrS %>%
                              dplyr::mutate(flrrm = format(LRR_mean_stab, digits = 2)) %>%
                              tidyr::unite(lrrstr, relation, flrrm, sep = ":"))$lrrstr, collapse  = " "))
  lrrMedians <- paste("LRR median",
                      paste0((lrrS %>%
                                dplyr::mutate(flrrmed = format(LRR_median_stab, digits = 1)) %>%
                                tidyr::unite(lrrstr, relation, flrrmed, sep = ":"))$lrrstr, collapse  = " "))
  lrrSds <- paste("LRR SD loc",
                  paste0((lrrS %>%
                            dplyr::mutate(flrrsd = format(LRR_SD_stab, digits = 2)) %>%
                            tidyr::unite(lrrstr, relation, flrrsd, sep = ":"))$lrrstr, collapse  = " "))

  lrrFSds <- paste("LRR SD flank",
                   paste0((lrrS %>%
                             dplyr::mutate(flrrsd = format(flank_LRR_SD_stab, digits = 2)) %>%
                             tidyr::unite(lrrstr, relation, flrrsd, sep = ":"))$lrrstr, collapse  = " "))

  # unite all the information
  ddt <- dbaf %>%
    dplyr::left_join(lrrS, by = c("relation")) %>%
    dplyr::left_join(flanktest, by = c("relation")) %>%
    dplyr::mutate(baftest = clbafset,
           lrrtest = cllrrset) %>%
    dplyr::mutate(sample = unique(input_data$sample),
           coordcnv = unique(input_data$coordcnv),
           type= unique(input_data$type),
           testLRR = lrrStr,
           distLRR = lrrDStr,
           meansLRR = lrrMeans,
           sdLRR = lrrSds,
           fsdLRR = lrrFSds,
           lrrmedtest = ifelse(cohen_cutoff==T,
                               ifelse(relation!="O",
                                      ifelse(lrr_wx_mu>sign_correct_baf, "cnvmd_at_zero",
                                             ifelse(lrr_wx_less<lrr_wx_gt, "cnv_less_flank", "cnv_gt_flank")),
                                      ifelse(lrr_wx_less<lrr_wx_gt, "cnv_less_flank", "cnv_gt_flank")),
                               "cnv_sim_flank"),
           lrrmedtest = ifelse(lrrmedtest == "cnv_sim_flank",
                               ifelse(lrr_wx_mu < 0.001, "cnv_mu_shift", lrrmedtest),
                               lrrmedtest),
           BAFprobes = ifelse(lrrmedtest == "cnv_gt_flank" & BAFprobes == "possiblyLOH",
                              ifelse(nprobes_50==0 & nprobes_dup >0, "BAFsuggDup", BAFprobes), BAFprobes)) %>%
    tidyr::unite(cnvTypeBAF, relation, BAFprobes, sep = ":", remove = FALSE) %>%
    dplyr::mutate(bafTypeTest = ifelse(type == "dup",
                                        ifelse(grepl("agreeNormal", BAFprobes[relation == "O"]), "BAFConflictsType", "BAFagreesType"),
                                        ifelse(grepl("possiblyLOH", BAFprobes[relation == "O"]), "BAFagreesType", "BAFConflictsType"))) %>%
    dplyr::mutate(perc_neg = N_NEG_LRR/nprobes*100,
                  perc_pos = N_POS_LRR/nprobes*100,
                  lrrTypeTest = ifelse(type == "dup",
                                       ifelse(perc_neg[relation == "O"] >=percent_opposite, "HighPercOppositeType", "LRRagreesType"),
                                       ifelse(perc_pos[relation == "O"] >=percent_opposite, "HighPercOppositeType", "LRRagreesType")))
  ddt <- ddt %>%
    dplyr::left_join(input_data %>%
                       dplyr::mutate(parental_orig=gsub("CN=[0-9] ", "", Support)) %>%
                       dplyr::filter(locus == "CNV") %>%
                       dplyr::select(sample, coordcnv, status, parental_orig, copynumber, numsnp) %>%
                       dplyr::distinct(), by = c("sample", "coordcnv")) %>%
    tidyr::separate(coordcnv, c("chr", "coord"), sep = ":", remove=F) %>%
    tidyr::separate(coord, c("start", "end")) %>%
    dplyr::mutate(length = as.integer(end) - as.integer(start))
  ddt
}
