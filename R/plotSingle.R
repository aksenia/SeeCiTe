#' Visualize local CNV SNP data and key summary statistics for a single individual
#'
#' @param input_data Data slot of an object read with readInputs, for a single CNV and individual
#' @param sifted_data Output of runAnalyzeSignal for the input data
#' @param print_title Print detailed summary statistics in a plot header. Default is TRUE.
#' @return A ggplot2 plot object
#' @export
#'
#' @examples
plotSingle <- function(input_data, sifted_data, print_title=TRUE) {
  ## Add global quality data if avaialble
  mono_marker <- T
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
    dplyr::mutate(color = ifelse(locus == "flank", locus, gsub("O", "?", relation)),
           color = factor(color, levels = c("?", "F", "M", "flank")))
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
    ggplot2::geom_point(ggplot2::aes(shape=locus)) +
    ggplot2::coord_cartesian(ylim = c(-limv, limv)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    myColIntensity[["color"]] +
    ggplot2::scale_shape_manual(values=c(19, 21) , guide = "none") +
    ggplot2::scale_alpha_manual(values=c(0.6, 0.8,0.8,0.9) , guide = "none") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                   axis.title.x = element_blank(),
          legend.position = "none") +
    ggplot2::ylab("Log R Ratio")

  gbaf <- ggplot2::ggplot(input_data %>% dplyr::filter(parameter == "BAF"),
                 ggplot2::aes(x = Position, y = value, color = color)) +
    ggplot2::geom_point(ggplot2::aes(shape=locus)) +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted") +
      ggplot2::geom_hline(yintercept = 0.33, linetype = "dotted") +
      ggplot2::geom_hline(yintercept = 0.66, linetype = "dotted") +
      ggplot2::scale_x_continuous(labels = scales::comma) +
    myColIntensity[["color"]] +
      ggplot2::scale_shape_manual(values=c(19, 21) , guide = "none") +
      ggplot2::scale_alpha_manual(values=c(0.6, 0.8,0.8,0.9) , guide = "none") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                     axis.title.x = element_blank(),
          legend.position = "none") +
      ggplot2:: ylab("B Allele Frequency")

  if(mono_marker) {
    gbaf <- gbaf +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::ylab("B Allele Frequency w/o mono-probes")
  }

  lrr_belowZero <- paste("N below 0",
                         paste0((sifted_data %>%
                                   dplyr::summarise(N_NEG_LRR = unique(N_NEG_LRR)) %>%
                                   dplyr::mutate(flrr = format(N_NEG_LRR, digits = 2)) %>%
                                   tidyr::unite(lrrstr, flrr, sep = ":"))$lrrstr, collapse  = " "))

  lrr_aboveZero <- paste("N above 0",
                         paste0((sifted_data %>%
                                   dplyr::summarise(N_POS_LRR = unique(N_POS_LRR)) %>%
                                   dplyr::mutate(flrr = format(N_POS_LRR, digits = 2)) %>%
                                   tidyr::unite(lrrstr, flrr, sep = ":"))$lrrstr, collapse  = " "))

  lrrCountsStr <- paste(lrr_belowZero, lrr_aboveZero, sep = " -- ")
  #  wilcox test p-values on the plot

  consistency_str <- paste0(sifted_data %>% dplyr::ungroup() %>%
                              dplyr::select(cnvTypeBAF) %>%
                              dplyr::distinct() %>% unlist(use.names = F), collapse = ", ")

  baf_string <- paste0("BAF status: ", paste0(sifted_data %>%
                         dplyr::ungroup() %>%
                         dplyr::select(BAFprobes) %>%
                         dplyr::distinct() %>%
                         dplyr::rename(label=BAFprobes) %>%
                         unlist(use.names = F), collapse = ", "))

  flank_string <- paste("CNV and flank differ(KDE test)",
                        paste0(sifted_data %>%
                                 dplyr::ungroup() %>%
                                 dplyr::select(flanktest, pvalFlank, cohen_cutoff, cohen) %>%
                                 dplyr::mutate(decision_shift = ifelse(cohen_cutoff==T, format(pvalFlank, digits =2), paste0("Cohen's d:", cohen)),
                                               flanktest = gsub("O", "", flanktest)) %>%
                                 dplyr::distinct() %>%
                                 tidyr::unite(label, flanktest, decision_shift , sep = ":") %>%
                                 dplyr::group_by(cohen_cutoff) %>%
                                 dplyr::summarise(str = paste(label, collapse = ", ")) %>%
                                 tidyr::unite(lbl, cohen_cutoff, str, sep = " ") %>%
                                 unlist(use.names = F),
                               collapse = " | "))

  ## distribution shift function
  sft <- calcShift(input_data %>%
                     dplyr::filter(parameter == "LRR"), single = TRUE)
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

####
  cnvlen <- sifted_data %>%
    dplyr::mutate(len = round(length/1000)) %>%
    dplyr::select(len) %>%
    unlist(use.names=F)

  cnvmarkers <- input_data %>%
    dplyr::filter(locus == "CNV") %>%
    dplyr::group_by(sample, coordcnv) %>%
    dplyr::summarise(nmarkers = dplyr::n_distinct(Name)) %>%
    dplyr::ungroup()  %>%
    dplyr::select(nmarkers) %>%
    unlist(use.names = F)

  lrr_string <- paste("Copynumber ", unique(sifted_data$copynumber), " | ", unique(sifted_data$meansLRR), " | ", unique(sifted_data$fsdLRR), sep = "")

  title_top <- paste(
    paste(paste(unique(input_data$sample), unique(input_data$coordcnv), sep = ", "),
          "|| Length", paste(cnvlen, "Kb", sep = ""), "Markers:", cnvmarkers),
    lrr_string,
    flank_string,
    paste(baf_string, lrrCountsStr, sep = "  -- "),
    paste("--", LRRconsistency, "--", bafTypeTest),
    sep = "\n")

  title_anonymous <- paste(
    paste(paste("Individual", unique(input_data$sample),
                                    unique(input_data$coordcnv), sep = ", "),
          "|| Length", paste(cnvlen, "Kb", sep = ""), "Markers:", cnvmarkers),  sep = "\n")

  title_print <- ifelse(isTRUE(print_title), title_top,  title_anonymous)
  lay <- rbind(c(1, 1, 1, 3, 3, 3), c(1, 1, 1, 3, 3, 3),
               +       c(2, 2, 2, 3, 3, 3), c(2, 2, 2, 3, 3, 3))
  options(warn=-1)
  suppressMessages(gridExtra::grid.arrange(glrr, gbaf, shiftpl,
                       layout_matrix = lay,
                       top = title_print))
  options(warn=0)
}
