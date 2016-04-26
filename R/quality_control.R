##' Quality Control Plot
##'
##' Plot the quality of each base against the reference. Allow a quick graphical
##' representation of the base calling quality.
##' @param data the dataset.
##' @param mut the mutant to filter by.
##' @return a ggplot2 plot
##' @author Samuel Barreto
##' @import dplyr
##' @import ggplot2
##' @importFrom assertthat assert_that is.flag
##' @export

quality_control <- function(data, mut = c("ws", "sw", "w", "s"), clean=TRUE)
{
    assert_that(
        is.data.frame(data),
        is.vector(mut),
        is.flag(clean)
    )

    if (clean) {
        data <- data %>% filter(mutant %in% mut) %>% keep_clean_only()
    } else {
        data <- data %>% filter(mutant %in% mut)
    }

    ggplot(data, aes(x = refp, y = qual, color = qual )) +
        geom_point(aes(alpha = qual), size = 1) +
        geom_line(aes(group = name), alpha = 1/10, size = 0.1 ) +
        scale_color_continuous(high = "white", low = "black") +
        scale_alpha(range = c(0.5, 0.01)) +
        labs(x = "Position", y = "Qualité")
}
