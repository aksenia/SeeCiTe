#' Annotate original PennCNV trio data table with inheritance status using trio state decomposition
#'
#' @param penntriotable A data.frame read into R with readPennTrioOriginal
#'
#' @return A data.frame with status of inheritance, one of "denovo"/"inherited"/"ambiguous" for offspring and "call_in_parent" for parents
#' @export
#'
extractInheritance <- function(penntriotable){
  # extract and annotate parental calls
  penntrio_parents <- penntriotable %>%
    dplyr::filter(relation!="offspring") %>%
    dplyr::mutate(status = "call_in_parent")

  # extract offspring calls
  penntrio_offspring <- penntriotable %>%
    dplyr::filter(relation=="offspring") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(len_tstate = length(unique(strsplit(triostate, split="-")[[1]])))

    # analyze separately
  # ambiguous and unambiguos state
  ambstate <- penntrio_offspring %>%
    dplyr::filter(len_tstate > 1) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(status = paste(sort(unique(sapply(strsplit(triostate, split="-")[[1]], statusFromTriostate))), collapse="-"),
           status=ifelse(grepl("-", status), "ambiguous", status))

  unambstate <- penntrio_offspring %>%
    dplyr::filter(len_tstate == 1)

  unamb_states <- unambstate %>%
    dplyr::select(triostate) %>%
    dplyr::distinct() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(status = statusFromTriostate(triostate))

  unambstate <- unambstate %>%
    dplyr::left_join(unamb_states, by = c("triostate"))
  # merge back
  penntrio_offspring_st <- do.call("rbind", list(ambstate, unambstate)) %>%
    dplyr::select(-len_tstate) %>%
    dplyr::distinct()
  ## some events are listed several times with different triostate
  ## need to collapse those
  penntrio_st <- do.call("rbind", list(penntrio_offspring_st[,colnames(penntrio_parents)],
                                       penntrio_parents)) %>%
    dplyr::distinct()
  penntrio_st
}
