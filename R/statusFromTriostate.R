#' Extract inheritance status from penn trio string collapsing copy number strings to discrete "del", "dup", "normal" or ""copy-neutrLOH" following PennCNV encoding (Wang et al, 2007)
#'
#' @param penntriostate A string of three digits representing copynumber state for members of a trio/quad
#' @param offspring_order Can be used both on trios and quads, therefore the index of an offspring (1 or 2) needs to be specified. Defaut is 1. (a trio)
#'
#' @return A string with a status, one of "denovo"/"inherited"
#' @export
#'
#' @examples statusFromTriostate("323")
statusFromTriostate <- function(penntriostate, offspring_order = 1){
  penn_cn_state <- decodeHMMstates()
  idx <- 2 + offspring_order
  dt <- data.frame(relation=c("father", "mother", "offspring"),
                   state=c(substr(penntriostate, 1, 1),
                           substr(penntriostate, 2, 2),
                           substr(penntriostate, idx, idx)),
                   stringsAsFactors = F) %>%
    dplyr::mutate(state=paste0("state", state))

  result <- dt %>%
    dplyr::left_join(penn_cn_state, by = c("state")) %>%
    dplyr::select(-state) %>%
    tidyr::spread(relation, type) %>%
    dplyr::mutate(status = ifelse(offspring !=father && offspring != mother, "denovo", "inherited"))
  result$status
}
