#' Write the pdf files with the visualisations of local CNV SNP data for each trio in the data
#'
#' @param main_data A list holding input data with one obligatory table "data" and two optional "qcsum" and "merge"
#' output of readInputs
#' @param sifted_data A table with summary statistis per individual in trio, output of runAnalyzeSignal
#' @param classified_data A table with final classifications, offspring-centered. Output of classifyTrios
#' @param output_dir Folder in which to write the plot files
#' @param dataset Dataset string to use in file naming
#' @param mono_marker (TRUE/FALSE) If the array contains non-polymorphic markers with BAF>1, do not plot these. TRUE by default
#' @param subset_nprobes Only plot the CNVs containing that many or more probes
#' @param subset_length Only plot the CNVs of the given legth (in bp) or longer
#' @return
#'
#' @examples
plotCohort <- function(main_data, sifted_data, classified_data, output_dir, dataset, mono_marker=T, subset_nprobes=NULL, subset_length=NULL){
  ## create a folder for each file being visualized
  dir.create(output_dir, showWarnings = FALSE)
  data_calls <- main_data[["data"]]
  classified_data <- classified_data %>%
    tidyr::unite(prefix, siftt, param, sep=".", remove=F)
  if (!is.null(subset_nprobes)) {
    sifted_data <- sifted_data %>%
      dplyr::filter(numsnp>=subset_nprobes)
    data_calls <- data_calls %>%
      dplyr::filter(coordcnv %in% unique(sifted_data$coordcnv))
    classified_data <- classified_data %>%
      dplyr::filter(coordcnv %in% unique(sifted_data$coordcnv))
  }
  if (!is.null(subset_length)) {
    sifted_data <- sifted_data %>%
      dplyr::filter(length>=subset_length)
    data_calls <- data_calls %>%
      dplyr::filter(coordcnv %in% unique(sifted_data$coordcnv))
    classified_data <- classified_data %>%
      dplyr::filter(coordcnv %in% unique(sifted_data$coordcnv))
  }
  l_class <- with(classified_data , split(classified_data, prefix))

  #####  OUTPUT
  sapply(names(l_class), function(nm) {
    print(paste0("Plotting file ", nm))
    cnvlist <- unique(l_class[[nm]]$coordcnv)
    # sort by chromosome
    cnvlist <- names(sort(sapply(cnvlist, function(v) as.integer(gsub("chr|:.*", "", v)))))
    fname <- paste("LogRR_BalleleFreq_", dataset, "_", nm, ".pdf", sep = "")
    fname <- file.path(output_dir, fname)

    pdf(file = fname, width = 16)
    lapply(cnvlist, function(cnv){
      dfa <- data_calls %>%
        dplyr::filter(parameter %in% c("LRR", "BAF"), coordcnv == cnv) %>%
        dplyr::arrange(coordcnv, sample)
      # if samples share coordinates of a cnv
      # but have different patterns of parameters
      # we want them in different files!
      lapply(intersect(unique(dfa$sample), unique((l_class[[nm]] %>% dplyr::filter(coordcnv == cnv))$sample)), function(s){
        print(paste0("Plotting CNV ", cnv, " in ", s))
        sdfa <- dfa %>%
          dplyr::filter(sample == s)
        if (sum(grepl("qcsum", names(main_data)))>0) {
          print(paste0("QC summaries supplied "))
          sqcsum <- main_data[["qcsum"]] %>%
            dplyr::filter(sample == s)
        } else {
          sqcsum <- NULL
        }
        if (sum(grepl("merge", names(main_data)))>0) {
          print(paste0("Merging trace supplied "))
          smerge <- main_data[["merge"]] %>%
            dplyr::filter(sample == s, coordcnv==cnv)
        } else {
          smerge <- NULL
        }
        ssifted <- sifted_data %>%
          dplyr::filter(sample == s, coordcnv == cnv)
        plotRawTrio(input_data = sdfa, sifted_data = ssifted,
                      penn_qcsum = sqcsum,
                      mono_marker = mono_marker,
                      merge_trace=smerge)

      })
    })
    dev.off()
  })
}
