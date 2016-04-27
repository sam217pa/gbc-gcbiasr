##' Plot the mean gc content per sequence
##'
##' @importFrom assertthat is.string assert_that
##' @importFrom ggthemes scale_colour_solarized extended_range_breaks
##' @import dplyr
##' @import ggplot2
##' @export

plot_gc_content <- function(data, mut=NULL, plot_title=NULL)
{
    assert_that(
        is.data.frame(data),
        is.string(mut) | is.null(mut),
        is.string(plot_title) | is.null(plot_title)
    )

                                        # set default parameter
    if (is.null(mut)) mut <- c("ws", "sw", "w", "s")
    if (is.null(plot_title)) plot_title <- "Comparaison des taux de GC"


    data <- keep_clean_only() %>%
        get_gc_content(mut = mut) %>%
        mutate(mutant = factor(mutant, levels = c("w", "s", "sw", "ws"),
                               labels = c("AT", "GC", "AT/GC", "GC/AT")))

    ggplot(data, aes(x = type, y = GC)) +
        geom_jitter(aes(color = mutant), width = 1/4, alpha = 1/3, size = 1) +
        geom_line(aes(group = name, color = mutant), alpha = 1/20) +
        stat_summary(
            fun.y = mean, geom = "point", color = "black", shape = "¦", size = 6) +
        facet_grid(mutant ~ .) +
        labs(y = "% GC", x = "", title = plot_title) +
        scale_color_solarized() +
        coord_flip()
}
