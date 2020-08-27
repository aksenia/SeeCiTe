#' Annotate original PennCNV trio data table with inheritance status using trio state decomposition
#'
#' @param filename_orig A path to output file of PennCNV trio
#' @param filename_merged A path to a result of merging with PennCNV of the filename_orig
#'
#' @return  Returns paths to and writes files with status of inheritance, one of "denovo"/"inherited"/"ambiguous" for offspring and "call_in_parent" for parents (with and without parental data),
#' as well as file with no inheritance information for further merging. In case if no merging is planned, a log file with original PennCNV trio states is also written for further in the pipeline.
#' @export
#' @examples
runExtractInheritance <- function(filename_orig, filename_merged){
  orig_dt <- readPennTrioOriginal(filename = filename_orig)
  merged_dt <- readPennTrioOriginal(filename = filename_merged, triostate = FALSE)

  merge_log <- mapTriostate(orig_dt = orig_dt,  merged_dt =  merged_dt)

  merged_dt <- merged_dt %>%
    dplyr::left_join(merge_log %>% dplyr::select(sample, coordcnv, triostate)) %>%
    dplyr::left_join(orig_dt %>% dplyr::select(sample, relation)) %>%
    dplyr::distinct()

  penntrio_status <- extractInheritance(penntriotable = merged_dt)
  out_file <- penntrio_status %>%
    dplyr::mutate(numsnp = paste("numsnp", numsnp, sep = "="),
           length = paste("length", length, sep = "="),
           cn = paste("cn", copynumber, sep = "="),
           startsnp = paste("startsnp", startsnp, sep = "="),
           endsnp = paste("endsnp", endsnp, sep = "="),
           triostate = paste("triostate", status, sep = "="),
           HMMstate=gsub("state1", "state2", state),
           HMMstate=gsub("state6", "state5", state),
           cn=gsub("0","1",cn),
           cn=gsub("4","3",cn)) %>%
    dplyr::select(-state) %>%
    tidyr::unite(state, HMMstate, cn, sep = ",")


  unmasked_status <- out_file %>%
    dplyr::select(coordcnv, numsnp, length, state, sample, startsnp, endsnp, relation, status)

  annotated_status <- out_file %>%
    tidyr::unite(sampleid, sample, relation, status, sep=":") %>%
    dplyr::select(coordcnv, numsnp, length, state, sampleid, startsnp, endsnp)

  ## define output
  prefix <- tools::file_path_sans_ext(filename_merged)
  out_unmasked <- paste(prefix, "_unmaskedstatus.triocnv", sep = "")
  out_annotated <- paste(prefix, "_annot.triocnv", sep = "")
  out_annotated_offspring <- paste(prefix, "_annot_offspring.triocnv", sep = "")
  out_mergelog <- paste(prefix, "_merge.log", sep = "")
  print(paste0("Writing files with prefix ", prefix))
  write.table(unmasked_status, file = out_unmasked, quote = F, col.names = F, row.names = F, sep = "\t")
  write.table(annotated_status, file = out_annotated, quote = F, col.names = F, row.names = F, sep = "\t")
  write.table(annotated_status %>% dplyr::filter(grepl("offspring", sampleid)), file = out_annotated_offspring, quote = F, col.names = F, row.names = F, sep = "\t")
  write.table(merge_log %>% dplyr::select(-triostate), file = out_mergelog, quote = F, col.names = F, row.names = F, sep = "\t")
  print("Done!")
  return(list(triocnv_file=out_annotated_offspring, merge_trace=out_mergelog))
}
