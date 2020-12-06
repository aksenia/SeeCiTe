#' Write the CNV files with SeeCiTe assessment and all summary metrics in txt, plink and bed formats
#'
#' @param classified_data classified_data A table with final classifications, offspring-centered. Output of classifyTrios.
#' @param output_dir output_dir Folder in which to write the files.
#' @param dataset dataset Dataset string to use in file naming.
#' @param subset_nprobes Only output the CNVs containing that many or more probes
#' @param subset_length Only output the CNVs of the given legth (in bp) or longer

#' @return No value returned. Writes files into a given folder.

#' @export
#' @examples
writeSeecite <- function(classified_data, output_dir, dataset, subset_nprobes=NULL, subset_length=NULL){
  dir.create(output_dir, showWarnings = FALSE)
  if (!is.null(subset_nprobes)) {
    classified_data <- classified_data %>%
      dplyr::filter(numsnp >= subset_nprobes)
  }
  if (!is.null(subset_length)) {
    classified_data <- classified_data %>%
      dplyr::filter(length >= subset_length)
  }
  classified_data <- classified_data %>%
    tidyr::unite(prefix, seecite, param, sep=".", remove=F)
  l_class <- with(classified_data , split(classified_data, prefix))


  sapply(names(l_class), function(nm) {
    cnvTable <- l_class[[nm]] %>%
      tidyr::separate(coordcnv, c("chr", "coord"), sep = ":", remove = F) %>%
      dplyr::mutate(chr = as.integer(gsub("chr", "", chr))) %>%
      dplyr::arrange(chr, coordcnv, sample) %>%
      dplyr::select(-dplyr::matches("diff_"), -param, -dplyr::matches("lrrm_dev"),
                    -flanktest, - offspring_flank_lrrsd_pass, -dplyr::matches("mosaic"),
                    -bafTypeDisagr, -countBafTypeDisagrP, -chr, -coord)
    print(paste0("Writing files ", nm, ", with ", nrow(cnvTable), " CNVs"))
    #### write summary table file
    fn <- paste("List_", dataset, ".", nm, ".txt", sep = "")
    fn <- file.path(output_dir, fn)
    write.table(cnvTable, file = fn, quote = F, row.names = F, col.names = T, sep = "\t")
    ### write plink files
    cnvPlink <- l_class[[nm]] %>%
      tidyr::separate(coordcnv, c("s", "BP2"), sep = "-", remove=F) %>%
      dplyr::mutate(BP1=gsub(".*:", "", s), FID=sample, SCORE=0, IID = paste(sample, "o", sep=":")) %>%
      dplyr::rename(TYPE=copynumber, CHR=Chr, SITES=numsnp) %>%
      dplyr::select(FID, IID, CHR, BP1, BP2, TYPE, SCORE, SITES)

    fnp <- paste("Plink_", dataset, ".", nm, ".txt", sep = "")
    fnp <- file.path(output_dir, fnp)
    write.table(cnvPlink, file = fnp, quote = F, row.names = F, col.names = T, sep = "\t")
    ### write bed files
    cnvBed <- l_class[[nm]] %>%
      tidyr::separate(coordcnv, c("s", "BP2"), sep = "-", remove=F) %>%
      dplyr::mutate(BP1=gsub(".*:", "", s),
             strand=".", dummy1="0", dummy2="0", rgb=ifelse(copynumber==1, "255,0,0", "0,255,0")) %>%
      dplyr::select(Chr, BP1, BP2, sample, numsnp, strand, dummy1, dummy2, rgb)

    fnb <- paste("Bed_", dataset, ".", nm, ".bed", sep = "")
    fnb <- file.path(output_dir, fnb)
    write.table(cnvBed, file = fnb, quote = F, row.names = F, col.names = F, sep = "\t")
  })
}
