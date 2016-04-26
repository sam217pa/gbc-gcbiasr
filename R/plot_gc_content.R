##' Plot the mean gc content per sequence
##'
##' @importFrom assertthat is.string assert_that
##' @importFrom seqinr GC
##' @importFrom tidyr gather
##' @importFrom ggthemes scale_colour_solarized extended_range_breaks
##' @import dplyr
##' @import ggplot2
##' @export

plot_gc_content <- function(data, mut=NULL, plot_title=NULL) {

    assert_that(
        is.data.frame(data),
        is.string(mut) | is.null(mut),
        is.string(plot_title) | is.null(plot_title)
    )

    if (is.null(mut)) mut <- c("ws", "sw", "w", "s")
    if (is.null(plot_title)) plot_title <- "Comparaison des taux de GC"

    get_gc_content <- function(data, mut) {
        ## clean_seq <- function(x) toString(x) %>% gsub(", ", "", . )

        gc_content <- data %>%
            mutate(snpb = as.character(snpb),
                   expb = as.character(expb),
                   refb = as.character(refb)) %>%
            filter(mutant %in% mut) %>%
            group_by(name, mutant) %>%
            ## compute gc content of the respective sequence
            summarise(exp = seqinr::GC(expb),
                      snp = seqinr::GC(snpb),
                      ref = seqinr::GC(refb))

        ## convert the data set to long form.
        gc_content %>%
            gather("type", "GC", 3:5) %>%
            mutate(type = factor(
                       type, levels = c("snp", "exp", "ref"),
                       labels = c("Donneur", "Recombinant", "Receveur")))
    }

    data <-
        filter(data, !(is.na(refb)), !(is.na(snpb)), !(is.na(expb))) %>%
        keep_clean_only() %>%
        get_gc_content(mut = mut)

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
