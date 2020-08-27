#' Visualize local CNV SNP data and key summary statistics for a trio
#'
#' @param input_data Data slot of an object read with readInputs, for a single CNV and individual
#' @param sifted_data Output of runAnalyzeSignal for the input data
#' @param print_title Print detailed summary statistics in a plot header. Default is TRUE.
#' @param penn_qcsum  Table with PennCNV QC summaries stored in "qcsum" by readInputs if supplied. Will add global LRR_SD and WF values if present. Default is NULL.
#' @param mono_marker (TRUE/FALSE) If the array contains non-polymorphic markers with BAF>1, do not plot these. TRUE by default.
#' @param merge_trace Table in a custom format with original PennCNV trio status and merging statistics.
#' @return A ggplot2 plot object
#'
#' @examples
plotRawTrio <- function(input_data, sifted_data, print_title=TRUE, penn_qcsum=NULL, mono_marker=TRUE, merge_trace=NULL) {
  ## Add global quality data if avaialble
  if (!is.null(penn_qcsum)){
    qcs_string <- penn_qcsum %>%
      dplyr::mutate(qc_global = paste0("Global LRR SD:",format(LRR_SD, digits=2), " | WF:", format(WF, digits=2))) %>%
      dplyr::select(qc_global) %>%
      unlist(use.names=F)
    qcs_string <- paste0("|| ", qcs_string)
  } else {
    qcs_string <- ""
  }

  # color code
  createColorScaleIntensity <- function(name){
    vals <- c("#D2AF81FF", "#F05C3BFF", "#71D0F5FF", "#075149FF", "#71D0F5FF")
    names(vals) <- c("flank", "F", "O", "M", "?")
    colscale <- list()
    colscale[["color"]] <- ggplot2::scale_colour_manual(name = name, values = vals)
    colscale[["fill"]] <- ggplot2::scale_fill_manual(name = name, values = vals)
    colscale
  }
  myColIntensity <- createColorScaleIntensity(name="color")

  input_data <- input_data %>%
    dplyr::mutate(color = ifelse(locus == "flank", locus,
                                 ifelse(relation == "O" & type=="del", Origin,
                                                                 gsub("O", "?", relation))),
                  Support=ifelse(type=="del", Support, gsub(" .*", "", Support)),
           color = factor(color, levels = c("?", "F", "M", "flank")),
           relation = factor(relation, levels = c("O", "F", "M")))
  # set the symmetric LRR y-axis limits
  suppressWarnings(deciles_df <- input_data %>%
    dplyr::filter(parameter == "LRR") %>%
    dplyr::tbl_df() %>%
    tidyr::nest(-relation) %>%
    dplyr::mutate(Quantiles = purrr::map(data, ~ quantile(.$value, probs = seq(0, 1, 0.1))),
           Quantiles = purrr::map(Quantiles, ~ dplyr::bind_rows(.) %>% tidyr::gather())) %>%
    tidyr::unnest(Quantiles))

  limv <- min(c(unique((deciles_df %>% dplyr::filter(grepl("100", key), value==max(value)))$value), 2))
  # main LRR scatterplot
  glrr<- ggplot2::ggplot(input_data %>% dplyr::filter(parameter == "LRR"),
                ggplot2::aes(x = Position, y = value, color = color)) +
    ggplot2::geom_point(ggplot2::aes(shape=locus, alpha = color)) +
    ggplot2::coord_cartesian(ylim = c(-limv, limv)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    myColIntensity[["color"]] +
    ggplot2::scale_shape_manual(values=c(19, 21) , guide = "none") +
    ggplot2::scale_alpha_manual(values=c(0.6, 0.8,0.8,0.9) , guide = "none") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
          legend.position = "top") +
    ggplot2::ylab("Log R Ratio") +
    ggplot2::facet_wrap(~ relation, ncol = 1)

  gbaf <- ggplot2::ggplot(input_data %>% dplyr::filter(parameter == "BAF"),
                 ggplot2::aes(x = Position, y = value, color = color)) +
    ggplot2::geom_point(ggplot2::aes(shape=locus, alpha = color)) +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted") +
      ggplot2::geom_hline(yintercept = 0.33, linetype = "dotted") +
      ggplot2::geom_hline(yintercept = 0.66, linetype = "dotted") +
      ggplot2::scale_x_continuous(labels = scales::comma) +
    myColIntensity[["color"]] +
      ggplot2::scale_shape_manual(values=c(19, 21) , guide = "none") +
      ggplot2::scale_alpha_manual(values=c(0.6, 0.8,0.8,0.9) , guide = "none") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
          legend.position = "top") +
      ggplot2:: ylab("B Allele Frequency") +
      ggplot2::facet_wrap(~ relation , ncol = 1)

  if(mono_marker) {
    gbaf <- gbaf +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::ylab("B Allele Frequency w/o mono-probes")
  }

  ## pairwise test results
  # precalculated in sifted_data
  test_string <- unique(sifted_data$distLRR)

  # fetch the pair
  pr <- unique(sifted_data$lrrtest)

  lrr_string <- paste(test_string,  ", pair ", pr, " has least dist. LRR | ",
                      unique(sifted_data$meansLRR), " | ", unique(sifted_data$fsdLRR), sep = "")


  lrr_belowZero <- paste("N below 0",
                         paste0((sifted_data %>%
                                   dplyr::group_by(relation) %>%
                                   dplyr::summarise(N_NEG_LRR = unique(N_NEG_LRR)) %>%
                                   dplyr::mutate(flrr = format(N_NEG_LRR, digits = 2)) %>%
                                   tidyr::unite(lrrstr, relation, flrr, sep = ":"))$lrrstr, collapse  = " "))

  lrr_aboveZero <- paste("N above 0",
                         paste0((sifted_data %>%
                                   dplyr::group_by(relation) %>%
                                   dplyr::summarise(N_POS_LRR = unique(N_POS_LRR)) %>%
                                   dplyr::mutate(flrr = format(N_POS_LRR, digits = 2)) %>%
                                   tidyr::unite(lrrstr, relation, flrr, sep = ":"))$lrrstr, collapse  = " "))

  lrrCountsStr <- paste(lrr_belowZero, lrr_aboveZero, sep = " -- ")
  #  wilcox test p-values on the plot
  my_comparisons <- list( c("F:CNV", "M:CNV") , c("O:CNV", "F:CNV"), c("O:CNV", "M:CNV"))
  # boxplot for CNV/flank
  msdlrr <- ggpubr::ggviolin(input_data %>%
                       dplyr::filter(parameter == "LRR") %>%
                       tidyr::unite(rel, relation, locus, sep = ":", remove = F), x = "rel", y = "value",
                     fill = "relation", add = "boxplot", add.params = list(fill = "white")) +
    myColIntensity[["fill"]] +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired = F, label = "p.signif", size=2.5) +
    ggplot2::theme(legend.position = "none",
          axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle=90),
          plot.margin = ggplot2::unit(c(0.5, 1, 0, 1), "lines")) +
    ggplot2::ylab("Log R Ratio")

  consistency_str <- paste0(sifted_data %>% dplyr::ungroup() %>%
                              dplyr::select(cnvTypeBAF) %>%
                              dplyr::distinct() %>% unlist(use.names = F), collapse = ", ")

  baf_string <- paste0(sifted_data %>%
                         dplyr::ungroup() %>%
                         dplyr::select(relation, BAFprobes) %>%
                         dplyr::distinct() %>%
                         dplyr::group_by(BAFprobes) %>%
                         dplyr::summarise(rel = paste(relation, collapse = ",")) %>%
                         tidyr::unite(label, rel, BAFprobes, sep = ":") %>%
                         unlist(use.names = F), collapse = ", ")

  flank_string <- paste("CNV and flank differ(KDE test)",
                        paste0(sifted_data %>%
                                 dplyr::ungroup() %>%
                                 dplyr::select(flanktest, pvalFlank, cohen_cutoff, cohen) %>%
                                 dplyr::mutate(decision_shift = ifelse(cohen_cutoff==T, format(pvalFlank, digits =2), paste0("Cohen's d:", cohen))) %>%
                                 dplyr::distinct() %>%
                                 tidyr::unite(label, flanktest, decision_shift , sep = ":") %>%
                                 dplyr::group_by(cohen_cutoff) %>%
                                 dplyr::summarise(str = paste(label, collapse = ", ")) %>%
                                 tidyr::unite(lbl, cohen_cutoff, str, sep = " ") %>%
                                 unlist(use.names = F),
                               collapse = " | "))

  ## distribution shift function
  sft <- calcShift(input_data %>%
                     dplyr::filter(parameter == "LRR"))
  shiftpl <- sft$graph
  shiftd <- sft$data

  shifts <- shiftd[[1]] %>%
    dplyr::summarise(pos=sum(CNV>0), neg=sum(CNV<0), tot=length(CNV)) %>%
    dplyr::mutate(perc_pos = pos/tot*100, perc_neg = neg/tot*100)

  LRRconsistency <- sifted_data$lrrTypeTest

  bfp <- (sifted_data %>%
            dplyr::filter(relation == "O") %>%
            dplyr::select(BAFprobes) %>%
            dplyr::distinct())$BAFprobes

  bafTypeTest <- ifelse(unique(input_data$type) == "dup",
                        ifelse(grepl("agreeNormal", bfp), "BAFConflictsType", "BAFagreesType"),
                        ifelse(grepl("possiblyLOH", bfp), "BAFagreesType", "BAFConflictsType"))
  ##
  ptn <- sifted_data %>%
    dplyr::ungroup() %>%
    dplyr::select(relation, BAFprobes) %>%
    dplyr::distinct() %>%
    dplyr::group_by(BAFprobes) %>%
    dplyr::summarise(rel = paste(relation, collapse = ","), n_mem = n())
####
  # Individual similarity based on LRR and BAF profiles of CNV loci
  ## Pair of individuals similar based on LRR in CNV locus
  lrrset <- strsplit(pr,split="")[[1]]
  ## Each individual has a different BAF pattern
  if (sum(ptn$n_mem > 1) == 0) {
    bafset <- ptn$rel
    bafpb <- (ptn %>% filter(rel %in% lrrset))$BAFprobes
    msg <- ifelse(length(intersect(bafset, lrrset))>1,
                  ifelse(sum(grepl("possiblyLOH", bafpb))>0,
                               paste("LRR shows similarity for", pr, "with consistent BAF"),
                         paste("LRR shows similarity for", pr, "but BAF disagree")),
                               "Neither LRR not BAF tests are conclusive")
  } else if (sum(ptn$n_mem > 1) == 0) {
    ## All three individuals have the same BAF pattern
    bafset <- strsplit(ptn$rel, split=",")[[1]]
    bafpb=ptn$BAFprobes
    msg <- ifelse(length(intersect(bafset, lrrset))>1,
                  paste("LRR shows similarity for", pr, "with consistent BAF"),
                  "LRR test is not conclusive")
  } else {
    ## Two of the three individuals share BAF pattern
    bfset<- (ptn %>% dplyr::filter(n_mem > 1))$rel
    bafset <- strsplit(bfset, split = ",")[[1]]
    bafpb <- (ptn %>% dplyr::filter(n_mem > 1))$BAFprobes
    msg <- ifelse(length(intersect(bafset, lrrset))>1,
                  paste("LRR shows similarity for", pr, "with consistent BAF"),
                  ifelse(sum(grepl("possiblyLOH", bafpb))>0,
                         paste("LRR shows similarity for", pr, "with consistent BAF"),
                         paste("LRR (", pr, ") and BAF (", gsub(",","", bfset), ") tests disagree", sep="")))
  }
####
  cnvlen <- sifted_data %>%
    dplyr::mutate(len = round(length/1000)) %>%
    dplyr::select(len) %>%
    unlist(use.names=F)

  cnvmarkers <- input_data %>%
    dplyr::filter(locus == "CNV") %>%
    dplyr::summarise(nmarkers = n_distinct(Name)) %>%
    unlist(use.names = F)

  #### merging trace
  if (!is.null(merge_trace)){
    mtr_string <- paste0("Penn trio: ", merge_trace %>% dplyr::select(tstate), " Merged seg: ", merge_trace %>% dplyr::select(mlog))
  } else {
    mtr_string <- ""
  }

  title_top <- paste(
    paste(paste(unique(input_data$sample), unique(input_data$coordcnv), sep = ", "),
          "|| Length", paste(cnvlen, "Kb", sep = ""), "Markers:", cnvmarkers, qcs_string),
    gsub(", $", "", paste(unique(input_data$Support), paste("Status:", unique(input_data$status)), mtr_string, sep = ", ")),
    lrr_string,
    flank_string,
    paste(baf_string, lrrCountsStr, sep = "  -- "),
    paste(msg, "--", LRRconsistency, "--", bafTypeTest),
    sep = "\n")

  title_anonymous <- paste(
    paste(paste("Offspring", unique(input_data$coordcnv), sep = ", "),
          "|| Length", paste(cnvlen, "Kb", sep = ""), "Markers:", cnvmarkers),
    paste(unique(input_data$Support), paste("Status:", unique(input_data$status)), sep = ", "),
    sep = "\n")

  title_print <- ifelse(print_title==T, title_top,  title_anonymous)
  lay <- rbind(c(1, 1, 2, 2, 3, 3), c(1, 1, 2, 2, 3, 3),
                 c(1, 1, 2, 2, 5, 5), c(1, 1, 2, 2, 5, 5), c(1, 1, 2, 2, 5, 5))
  options(warn=-1)
  suppressMessages(gridExtra::grid.arrange(glrr, gbaf, msdlrr, shiftpl,
                       layout_matrix = lay,
                       top = title_print))
  options(warn=0)
}
