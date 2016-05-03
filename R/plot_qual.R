#' Plot the quality of the sequence
#'
#' @title plot_qual
#' @param data the output of \code{run_phruscle}.
#' @param mutant_ the mutant to filter by.
#' @return a ggplot2 plot.
#' @author Samuel Barreto
#' @importFrom RColorBrewer brewer.pal
#' @import dplyr
#' @import ggplot2
#' @export

plot_qual <- function(data, mutant_) {
    red    <- brewer.pal(n = 4, "Set1")[1]
    blue   <- brewer.pal(n = 4, "Set1")[2]
    grey   <- brewer.pal(n = 9, "Greys")[4]
    black  <- brewer.pal(n = 9, "Greys")[9]

    data %>%
        filter(mutant == mutant_) %>%
        ggplot(aes(x = refp, y = qual)) +
        geom_point(aes(color = qual), size = 1, alpha = 1/20) +
        geom_smooth(se = FALSE, color = red, linetype = "dotted") +
        ## scale_color_brewer(palette = "Greys") +
        scale_color_continuous(high = grey, low = black) +
        labs(x = "Position sur la référence", y = "")
}


#' Plot the alignment and the quality of the sequence juxtaposed
#'
#'
#' @title plot_align_qual
#' @param data the output of \code{run_phruscle}.
#' @param mutant_ the mutant to filter by.
#' @return a ggplot2 plot.
#' @author Samuel Barreto
#' @importFrom cowplot plot_grid
#' @export

plot_align_qual <- function(data, mutant_)
{
    plot_grid(
        plot_align(data, mutant_),
        plot_qual(data, mutant_),
        ncol = 1, rel_heights = c(8.5, 1.5),
        align = 'v', labels = c("  1", "  2")
    )
}
