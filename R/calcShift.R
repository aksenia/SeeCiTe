#' Create a decile plot and a scatter for CNV versus flanking region in an offspring. Uses rogme package.
#' @param lrr_dt A data table with probe level LRR values for all individuals in a trio.
#'
#' @return A list with a data table with sampling results and a graph object
#'
#' @import ggplot2
#'
#' @examples
calcShift <- function(lrr_dt){
  # we are interested on the offspring locus
  adf <- lrr_dt %>%
    dplyr::filter(relation == "O") %>%
    dplyr::ungroup() %>%
    dplyr::select(locus, value)  %>%
    dplyr::rename(obs = value, gr = locus) %>%
    dplyr::arrange(gr) %>%
    dplyr::mutate(gr = factor(gr, levels = c("CNV", "flank")))
  # for now just subset either of the set
  # subset the data.frame to same number of observations
  # in each category
  sadf <- do.call("rbind",
                  lapply(names(table(adf$gr)), function(n) {
                    t <- adf %>%
                      dplyr::filter(gr==n);
                    s = sample(1:nrow(t), min(table(adf$gr)));
                    t[s,]}))
  tsf <- rogme::shiftdhd(data = sadf,
                  formula = obs ~ gr, nboot = 200)
  # scatterplot
  p <- rogme::plot_scat2(adf,
                  xlabel = "",
                  ylabel = "",
                  alpha = .8,
                  shape = 21,
                  colour = "#71D0F5FF",
                  fill = "grey90") +
    ggplot2::geom_hline(yintercept = 0, linetype="dashed")

    p <- p + ggplot2::coord_flip() +
      ggplot2::ylab("CNV vs flanking loci in O") +
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = 10),
          axis.title.y = ggplot2::element_text(size = 10))
  # deciles
  psf <- rogme::plot_sf(tsf, plot_theme = 2, symb_size = 1.5)
  psf <- psf[[1]] +
    ggplot2::labs(x = paste(names(tsf)[1], "LRR deciles (O)"),
         y = paste("Diff deciles")) +
    ggplot2::theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  pl <- cowplot::plot_grid(p, psf, ncol = 1, nrow = 2,
                             rel_heights = c(1, 1), label_size = 10,
                             hjust = -0.5,
                             scale=.95)

  list(graph=pl, data=tsf)
}
