#' Link inheritance status and original trio state as well as merging tracing
#'
#' @param orig_dt A table with output file of PennCNV trio
#' @param merged_dt A Table with the result of merging with PennCNV of the original file
#'
#' @return A data table with merged sample and CNV ids as well as original tristates and counts of segments merged
#' @examples 

mapTriostate <- function(orig_dt, merged_dt){
  # read in the files
  # original data
  orig_list <- with(orig_dt, split(orig_dt, sample))
  orig_list <- lapply(orig_list, function(u) u$coordcnv)
  # merged data 
  merged_list <- with(merged_dt, split(merged_dt, sample))
  merged_list <- lapply(merged_list, function(u) u$coordcnv)
  out <- do.call("rbind", lapply(names(merged_list), function(n){
    a <- merged_list[[n]]
    b <- orig_list[[n]]
    a.sort <- bedr::bedr.sort.region(a);
    b.sort <- bedr::bedr.sort.region(b);
    bedr::bedr.join.region(a.sort, b.sort) %>%
    dplyr::mutate(coordcnv=paste0(V4,":", V5,"-", V6),
                  sample=n) %>%
    dplyr::select(index, coordcnv, sample)
    })) %>%
    dplyr::left_join(orig_dt) %>%
    dplyr::group_by(index, sample) %>%
    dplyr::summarise(tstate=paste0(sort(triostate), collapse="|"),
                     triostate=paste0(sort(triostate), collapse="-"),
                     mlog=n()) %>%
    dplyr::rename(coordcnv=index) %>%
    dplyr::arrange(sample, coordcnv)
  out
}






