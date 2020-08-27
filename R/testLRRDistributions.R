#' Apply a row of statistical tests to LRR distribtuions in CNV and flanking regions within and between individuals in a trio
#'
#' @param my_lrr The input data table with LRR values in CNV and flanking regions collated with readInputs in "data" for each individual
#' @param cohen_dist A cutoff on Cohen distance between CNV and flank LRR value distributions. If cohen distance is equal or smaller that this cutoff the distributions are considered not different without further testing. Default is 3.
#' @param sign_lrr_dist A cutoff of p-value at which difference between two individual's CNV LRR value distributions is set as significant, using Kernel density based global two-sample comparison test from ks package. Default is 0.01.
#'
#' @return A data table with summaries and test results for LRR-based tests between the individuals (lrrtest) and between CNV and flank for each individual (flanktest).
#'
#' @examples
testLRRDistributions <- function(my_lrr, cohen_dist=3, sign_lrr_dist=0.01){
  # compare LRR values distributions
  # Cohen's distance
  bcd <- with(my_lrr, list(FO=NA,
                       MO=NA,
                       FM=NA,
                       Oflank=effsize::cohen.d(d=na.omit(O), f=na.omit(flank_O), hedges.correction = T)$magnitude,
                       Fflank=effsize::cohen.d(d=na.omit(Fa), f=na.omit(flank_F), hedges.correction = T)$magnitude,
                       Mflank=effsize::cohen.d(d=na.omit(M), f=na.omit(flank_M), hedges.correction = T)$magnitude))

  # Kernel Density Smoothing
  bks<- with(my_lrr, list(FO=ks::kde.test(c(na.omit(Fa)), c(na.omit(O)))$pvalue,
                      MO=ks::kde.test(c(na.omit(M)), c(na.omit(O)))$pvalue,
                      FM=ks::kde.test(c(na.omit(Fa)), c(na.omit(M)))$pvalue,
                      Oflank=ks::kde.test(c(na.omit(O)), c(stabilizeDistr(flank_O)))$pvalue,
                      Fflank=ks::kde.test(c(na.omit(Fa)), c(stabilizeDistr(flank_F)))$pvalue,
                      Mflank=ks::kde.test(c(na.omit(M)), c(stabilizeDistr(flank_M)))$pvalue))

  # Kolmogorov-Smirnov Tests
  bks1<- with(my_lrr, list(FO=stats::ks.test(stabilizeDistr(O), stabilizeDistr(Fa), exact = FALSE, alternative = "two.sided")$p.value,
                       MO=stats::ks.test(stabilizeDistr(O), stabilizeDistr(M), exact = FALSE, alternative = "two.sided")$p.value,
                       FM=stats::ks.test(stabilizeDistr(Fa), stabilizeDistr(M), exact = FALSE, alternative = "two.sided")$p.value,
                       Oflank=stats::ks.test(stabilizeDistr(O), na.omit(flank_O), exact = FALSE, alternative = "two.sided")$p.value,
                       Fflank=stats::ks.test(stabilizeDistr(Fa), na.omit(flank_F), exact = FALSE, alternative = "two.sided")$p.value,
                       Mflank=stats::ks.test(stabilizeDistr(M), na.omit(flank_M), exact = FALSE, alternative = "two.sided")$p.value))
  # Mann-Whitney
  bmw<- with(my_lrr, list(FO=stats::wilcox.test(stabilizeDistr(O), stabilizeDistr(Fa), exact = FALSE, paired = FALSE)$p.value,
                      MO=stats::wilcox.test(stabilizeDistr(O), stabilizeDistr(M), exact = FALSE, paired = FALSE)$p.value,
                      FM=stats::wilcox.test(stabilizeDistr(Fa), stabilizeDistr(M), exact = FALSE, paired = FALSE)$p.value,
                      Oflank=stats::wilcox.test(stabilizeDistr(O), flank_O, exact = FALSE, paired = FALSE)$p.value,
                      Fflank=stats::wilcox.test(stabilizeDistr(Fa), flank_F, exact = FALSE, paired = FALSE)$p.value,
                      Mflank=stats::wilcox.test(stabilizeDistr(M), flank_M, exact = FALSE, paired = FALSE)$p.value))
  # Brunner Munzel (Generalized Wilcoxon)
  bbm <- with(my_lrr, list(FO= lawstat::brunner.munzel.test(stabilizeDistr(O), stabilizeDistr(Fa))$p.value,
                       MO= lawstat::brunner.munzel.test(stabilizeDistr(O), stabilizeDistr(M))$p.value,
                       FM=lawstat::brunner.munzel.test(stabilizeDistr(Fa), stabilizeDistr(M))$p.value,
                       Oflank=lawstat::brunner.munzel.test(stabilizeDistr(O), flank_O)$p.value,
                       Fflank=lawstat::brunner.munzel.test(stabilizeDistr(Fa), flank_F)$p.value,
                       Mflank=lawstat::brunner.munzel.test(stabilizeDistr(M), flank_M)$p.value))
  # Wilcox less
  bmwl<- with(my_lrr, list(Oflank=stats::wilcox.test(stabilizeDistr(O), flank_O, exact = FALSE, paired = FALSE, alternative = "less")$p.value,
                       Fflank=stats::wilcox.test(stabilizeDistr(Fa), flank_F, exact = FALSE, paired = FALSE, alternative = "less")$p.value,
                       Mflank=stats::wilcox.test(stabilizeDistr(M), flank_M, exact = FALSE, paired = FALSE, alternative = "less")$p.value))
  # Wilcox greater
  bmwg<- with(my_lrr, list(Oflank=stats::wilcox.test(stabilizeDistr(O), flank_O, exact = FALSE, paired = FALSE, alternative = "greater")$p.value,
                       Fflank=stats::wilcox.test(stabilizeDistr(Fa), flank_F, exact = FALSE, paired = FALSE, alternative = "greater")$p.value,
                       Mflank=stats::wilcox.test(stabilizeDistr(M), flank_M, exact = FALSE, paired = FALSE, alternative = "greater")$p.value))

  # Wilcox one-sample
  bmwm<- with(my_lrr, list(Oflank=stats::wilcox.test(stabilizeDistr(O), mu = 0, alternative = "two.sided")$p.value,
                       Fflank=stats::wilcox.test(stabilizeDistr(Fa), mu =0, alternative = "two.sided")$p.value,
                       Mflank=stats::wilcox.test(stabilizeDistr(M), mu =0, alternative = "two.sided")$p.value))

  bmwmg<- with(my_lrr, list(Oflank=stats::wilcox.test(O, mu = 0, alternative = "greater")$p.value,
                        Fflank=stats::wilcox.test(Fa, mu =0, alternative = "greater")$p.value,
                        Mflank=stats::wilcox.test(M, mu =0, alternative = "greater")$p.value))

  bmwml<- with(my_lrr, list(Oflank=stats::wilcox.test(O, mu = 0, alternative = "less")$p.value,
                        Fflank=stats::wilcox.test(Fa, mu =0, alternative = "less")$p.value,
                        Mflank=stats::wilcox.test(M, mu =0, alternative = "less")$p.value))


  hell <- sapply(list(c("O", "Fa"), c("O", "M"), c("Fa", "M"),
                      c("O", "flank_O"), c("Fa", "flank_F"), c("M", "flank_M")), function(n){
                        x <- na.omit(my_lrr[, n[[1]]]) %>% unlist(use.names=F)
                        x <- stabilizeDistr(x)
                        y <- na.omit(my_lrr[, n[[2]]]) %>% unlist(use.names=F)
                        if (!grepl("flank", n[[2]])){
                          y <- stabilizeDistr(y)
                        }
                        t <- try(statip::hellinger(x, y, method = 2))
                        if (class(t) == "try-error")
                          t <- try(statip::hellinger(x, y, method = 1))
                        t
                      }, simplify = F, USE.NAMES = T)


  result <- data.frame(ks_pval=unlist(bks),
                      cohen = unlist(bcd),
                      mw_pval=unlist(bmw),
                      bm_pval=unlist(bbm),
                      hellinger=unlist(hell),
                      lrrtest = names(bmw),
                      lrr_wx_less = c(rep(NA, 3), unlist(bmwl)),
                      lrr_wx_gt = c(rep(NA, 3), unlist(bmwg)),
                      lrr_wx_mu = c(rep(NA, 3), unlist(bmwm)),
                      lrr_wx_mu_gt = c(rep(NA, 3), unlist(bmwmg)),
                      lrr_wx_mu_less = c(rep(NA, 3), unlist(bmwml)),
                      stringsAsFactors = FALSE)
  # sometimes Br-M test gives NaN - detect it and replace by Wicx rank
  nancount <- with(result, sum(is.nan(bm_pval)))
  cd <- rlang::quo(cohen_dist)
  result <- result %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pvalLRR=ifelse(nancount==0, bm_pval, mw_pval),
           test=ifelse(nancount==0, "BrM.", "MaWh."),
           cohen_cutoff = cohen > cohen_dist,
           ks_cutoff = ks_pval < sign_lrr_dist) %>%
    dplyr::arrange(pvalLRR)
  result
}
