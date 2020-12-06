#' Apply a row of statistical tests to LRR distribuions in CNV and flanking regions within a single individual
#'
#' @param my_lrr The input data table with LRR values in CNV and flanking regions collated with readInputs in "data" for each individual
#' @param cohen_dist A cutoff on Cohen distance between CNV and flank LRR value distributions. If cohen distance is equal or smaller that this cutoff the distributions are considered not different without further testing. Default is 3.
#'
#' @return A data table with summaries and test results for LRR-based tests between the individuals (lrrtest) and between CNV and flank for each individual (flanktest).
#'
#' @examples
testLRRDistributionsSingle <- function(my_lrr, cohen_dist=3){
  # compare LRR values distributions
  # Cohen's distance
  bcd <- with(my_lrr, list(Oflank=as.numeric(effsize::cohen.d(d=na.omit(O), f=na.omit(flank_O), hedges.correction = T)$magnitude)))

  # Kernel Density Smoothing
  bks<- with(my_lrr, list(Oflank=ks::kde.test(c(na.omit(O)), c(stabilizeDistr(flank_O)))$pvalue))

  # Kolmogorov-Smirnov Tests
  bks1<- with(my_lrr, list(Oflank=stats::ks.test(stabilizeDistr(O), na.omit(flank_O), exact = FALSE, alternative = "two.sided")$p.value))
  # Mann-Whitney
  bmw<- with(my_lrr, list(Oflank=stats::wilcox.test(stabilizeDistr(O), flank_O, exact = FALSE, paired = FALSE)$p.value))
  # Brunner Munzel (Generalized Wilcoxon)
  bbm <- with(my_lrr, list(Oflank=lawstat::brunner.munzel.test(stabilizeDistr(O), flank_O)$p.value))
  # Wilcox less
  bmwl<- with(my_lrr, list(Oflank=stats::wilcox.test(stabilizeDistr(O), flank_O, exact = FALSE, paired = FALSE, alternative = "less")$p.value))
  # Wilcox greater
  bmwg<- with(my_lrr, list(Oflank=stats::wilcox.test(stabilizeDistr(O), flank_O, exact = FALSE, paired = FALSE, alternative = "greater")$p.value))

  # Wilcox one-sample
  bmwm<- with(my_lrr, list(Oflank=stats::wilcox.test(stabilizeDistr(O), mu = 0, alternative = "two.sided")$p.value))

  bmwmg<- with(my_lrr, list(Oflank=stats::wilcox.test(O, mu = 0, alternative = "greater")$p.value))

  bmwml<- with(my_lrr, list(Oflank=stats::wilcox.test(O, mu = 0, alternative = "less")$p.value))

  result <- data.frame(ks_pval=unlist(bks),
                      cohen = unlist(bcd),
                      mw_pval=unlist(bmw),
                      bm_pval=unlist(bbm),
                      lrrtest = names(bmw),
                      lrr_wx_less = unlist(bmwl),
                      lrr_wx_gt = unlist(bmwg),
                      lrr_wx_mu = unlist(bmwm),
                      lrr_wx_mu_gt = unlist(bmwmg),
                      lrr_wx_mu_less = unlist(bmwml),
                      stringsAsFactors = FALSE)
  # sometimes Br-M test gives NaN - detect it and replace by Wicx rank
  nancount <- with(result, sum(is.nan(bm_pval)))
  result <- result %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pvalLRR=ifelse(nancount==0, bm_pval, mw_pval),
           test=ifelse(nancount==0, "BrM.", "MaWh."),
           cohen_cutoff = cohen > cohen_dist,
           ks_cutoff = NA) %>%
    dplyr::arrange(pvalLRR)
  result
}
