#' Create a table mapping PennCNV HMM states to copy numbers
#'
#' @return A data.frame with the mapping of HMM states in PennCNV to copynumber states
#' @export
#'
#' @examples decodeHMMstates()
decodeHMMstates <- function(){
  data.frame(state = paste("state", c(1:6), sep = ""),
                            type = c(rep("del", 2), "normal", "copy-neutrLOH", rep("dup", 2)),
                            stringsAsFactors = F)
}
