#' Cluster BAF values using given cluster centers, helper function
#'
#' @param mybaf The input data table collated with readSiftInput for individual CNV and sample in "data".
#' @param margin.dup A margin value around dup clusters AAB or ABB (0.33 and 0.66). Default is 0.05.
#' @param margin.het A margin value around heterozygous AB cluster (0.5). Default is 0.1.
#' @param margin.nh A margin value around homozygous AA or BB clusters (0 and 1). Default is 0.17
#' @return A data table with counts and fractions of probes in each cluster.
#'
#' @examples
clusterBAF <- function(mybaf, margin.dup=0.05, margin.het=0.1, margin.nh=0.17){

  checkDupClusters <- function(cdf, margin = margin.dup){
    cdf %>%
      dplyr::ungroup() %>%
      dplyr::filter(!value %in% outliers::outlier(value)) %>%
      dplyr::summarise(nprobes_33 = sum(value > 0.33 - margin & value < 0.33 + margin),
                       nprobes_66 = sum(value > 0.66 - margin & value < 0.66 + margin),
                       nprobes_cn4 = sum((value > 0.25 -margin & value < 0.25 + margin)|(value > 0.75 -margin & value < 0.75 +margin)),
                       # the middle centers for cn5 are not informative as they overlap with normal state cluster
                       # count only the characteristic outliers on the sides
                       nprobes_cn5 = sum((value > 0.20 -margin & value < 0.20 + margin)|(value > 0.80 -margin & value < 0.80 +margin))) %>%
      unlist()
  }

  checkHetCluster <- function(cdf, margin = margin.het){
    cdf %>%
      dplyr::ungroup() %>%
      dplyr::summarise(nprobes_50 = sum(value > 0.50 - margin & value < 0.50 + margin)) %>%
      unlist()
  }

  checkNonHomCluster <- function(cdf, margin = margin.nh){
    cdf %>%
      dplyr::ungroup() %>%
      dplyr::summarise(nprobes_middle = sum(value >= 0 + margin & value <= 1 - margin)) %>%
      unlist()
  }

  result <- do.call("rbind", lapply(unique(mybaf$relation),  function(cn) {
    as <- mybaf %>%
      dplyr::filter(!is.na(value), relation == cn)

    dupClcount <- checkDupClusters(as)
    hetClcount <- checkHetCluster(as)
    hetClcount_precise <- checkHetCluster(as, margin = 0.04)
    midClcount <-  checkNonHomCluster(as)

    as <- as %>%
      dplyr::summarise(relation = cn,
             nprobes33_66 = dupClcount[["nprobes_33"]] + dupClcount[["nprobes_66"]],
             nprobes_cn4 = dupClcount[["nprobes_cn4"]],
             nprobes_cn5 = dupClcount[["nprobes_cn5"]],
             nprobes_50 = hetClcount[["nprobes_50"]],
             nprobes_50_center = hetClcount_precise[["nprobes_50"]],
             nprobes_dup = max(nprobes33_66, nprobes_cn4, nprobes_cn5),
             nprobes_middle = midClcount[["nprobes_middle"]],
             nprobes_dup_total = nprobes_middle - nprobes_50,
             nprobes_mono = length(value>1),
             nprobes = dplyr::n_distinct(Name)) %>%
      dplyr::ungroup()
    as
  }))
  ## we want to divide by only middle probes
  ## because there can be many hom probes with dup cluster in the middle!
  result %>%
    dplyr::group_by(relation) %>%
    dplyr::mutate(frac_50 = ifelse(nprobes_middle==0, 0, nprobes_50/nprobes_middle),
           frac_33_66 = ifelse(nprobes_middle==0, 0, nprobes33_66/nprobes_middle),
           frac_dup = ifelse(nprobes_middle==0, nprobes_dup_total, nprobes_dup_total/nprobes_middle),
           frac_mono = nprobes_mono/nprobes) %>%
    dplyr::ungroup()
}
