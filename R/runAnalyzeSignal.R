#' A wrapper around analyzeSignalConcordance function to run for a cohort of trios
#'
#' @param input_data The input data table ("data") collated with readInputs.
#' @param args Global arguments for the dataset analysis run.
#' @param use_cache Use pre-calculated object from previous runs which is written every time the function is run into the folder defined by "cache_id" in args. This analysis takes a while on a large dataset. If it is necessary to re-run from scratch set this argument to FALSE. Default is FALSE.
#' @return A data table with SeeCiTe preparation analysis for all trios and CNVs in the input.
#' @export
#' @examples
runAnalyzeSignal <- function(input_data, args, use_cache=FALSE){
  cache_fn <- file.path(args[["cache_id"]], paste(args[["dataset"]], "_clu_baf", ".rds", sep = ""))
  if (use_cache ==TRUE && file.exists(cache_fn)) {
    print(paste0("Found cached file ", cache_fn, ", reading in"))
    result <- readRDS(file = cache_fn)
  } else {
    print(paste0("No cache supplied , doing the analysis from scratch"))
    result <- do.call("rbind", lapply((input_data %>%
                      dplyr::arrange(Chr, Position) %>%
                      dplyr::select(coordcnv) %>%
                      unique() %>%
                      unlist(use.names = F)), function(cnv) {
                        print(cnv)
                        dfa <- input_data %>%
                          dplyr::filter(parameter %in% c("LRR", "BAF"), coordcnv == cnv) %>%
                          dplyr::mutate(color = ifelse(relation == "O", Origin, relation))
                        ss <- lapply(unique(dfa$sample), function(s){
                          print(s)
                          ## original data
                          sdfa <- dfa %>%
                            dplyr::filter(sample == s)
                            asdfa <- analyzeSignalConcordance(sdfa)
                            })
                        r <- do.call("rbind", ss)
                   }))
  saveRDS(result, file=cache_fn)
  }
  result
}

