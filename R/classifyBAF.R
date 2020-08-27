#' Cluster BAF values using given cluster centers, a helper function
#'
#' @param mybaf The input data table collated with readSiftInput for individual CNV and sample in "data"
#' @param margin.dup A margin value around dup clusters AAB or ABB (0.33 and 0.66). Default is 0.05.
#' @param margin.het A margin value around heterozygous AB cluster (0.5). Default is 0.1.
#' @param margin.nh A margin value around homozygous AA or BB clusters (0 and 1). Default is 0.17
#' @return A data table with counts and fractions of probes in each cluster defined by the arguments and copynumber status derived based on those, one of "BAFagreeNormal", "possiblyLOH" or "BAFsuggDup".
#'
#' @examples
classifyBAF <- function(mybaf, margin.dup=0.05, margin.het=0.1, margin.nh=0.17){
  ## for non-polymorphic probes (eg cnv...  or CN for Affy)
  ## with BAF=2, at the moment, exclude those
  ## need to handle a case when all probes are non-polymorphic

  dbaf <- clusterBAF(mybaf, margin.dup=margin.dup, margin.het=margin.het, margin.nh=margin.nh) %>%
    dplyr::mutate(norm_probes = nprobes_50 > 0,
                  dup_probes = nprobes33_66 > 0,
                  middle_probes = nprobes_middle/nprobes) %>%
    dplyr::mutate(BAFprobes = ifelse((middle_probes < 0.05 & nprobes_50 + nprobes33_66 <=2) | (nprobes_middle == 1) |
                                       (nprobes_50 + nprobes_dup<=1) |
                                       (nprobes_50==0 & nprobes_dup == 1 & nprobes_cn5==1), "possiblyLOH",
                                     ifelse(nprobes_50_center >= nprobes_dup_total,
                                            ifelse((nprobes_50_center ==  nprobes_dup_total & nprobes_dup==nprobes_cn4 & nprobes_dup > nprobes33_66 & nprobes_dup > 1 & nprobes_dup >nprobes_cn5) |
                                                     (nprobes_50_center ==  nprobes_dup_total & nprobes_dup == nprobes_cn5 & frac_dup > 0.25 & nprobes_dup > 1 & nprobes_50_center < nprobes_cn5),
                                                   "BAFsuggDup", "BAFagreeNormal"),
                                            ifelse(frac_dup < 0.1 | (nprobes_50 > nprobes_dup_total), "BAFagreeNormal", "BAFsuggDup"))))
  dbaf
}
