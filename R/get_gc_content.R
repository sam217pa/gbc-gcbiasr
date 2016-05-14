
#' compute the gc content per sequence
#'
#' @param data a dataframe output of \code{read_phruscle}
#' @param mut the mutant to filter by.
#' @param clean must clean the data before. Default is TRUE
#' @importFrom seqinr GC
#' @importFrom assertthat assert_that is.string
#' @importFrom tidyr gather
#' @import dplyr
#'
#' @export

get_gc_content <- function(data, mut=NULL, clean = TRUE) {

    assert_that(
        is.data.frame(data),
        is.null(mut) | is.vector(mut) | is.string(mut),
        is.flag(clean)
    )

    if (is.null(mut)) mut <- c("ws", "sw", "w", "s")
    if (clean) gc_content <- keep_clean_only(data)
    else gc_content <- data

    gc_content <- gc_content %>%
        ## omit the rare cases of N in sequences
        filter(!(is.na(refb)), !(is.na(snpb)), !(is.na(expb))) %>%
        mutate(snpb = as.character(snpb),
               expb = as.character(expb),
               refb = as.character(refb)) %>%
        filter(mutant %in% mut) %>%
        group_by(name, mutant) %>%
        ## compute gc content of the respective sequence
        summarise(exp = seqinr::GC(expb),
                  snp = seqinr::GC(snpb),
                  ref = seqinr::GC(refb)) %>%
        ## convert the data set to long form.
        gather("type", "GC", 3:5) %>%
        mutate(type = factor(
                   type, levels = c("snp", "exp", "ref"),
                   labels = c("Donneur", "Recombinant", "Receveur")))

    class(gc_content) <- c("gc_content", class(gc_content))
    gc_content
}

##' Plot the mean gc content per sequence
##'
##' @import dplyr
##' @import ggplot2
##' @export

plot.gc_content <- function(x, plot_title=NULL, ...)
{
    ## set default parameter
    if (is.null(plot_title)) plot_title <- "Comparaison des taux de GC"

    x %>%
        mutate(
            mutant = factor(mutant, levels = c("w", "s", "sw", "ws"),
                            labels = c("AT", "GC", "AT/GC", "GC/AT"))
        ) %>%
    ggplot(aes(x = type, y = GC, color = type)) +
        geom_jitter(width = 1/4, alpha = 1/3, size = 1) +
        geom_line(aes(group = name), alpha = 1/20) +
        stat_summary(
            fun.y = mean
           ,geom = "point", color = "black", shape = "¦", size = 6) +
        facet_grid(mutant ~ .) +
        labs(y = "% GC", x = "", title = plot_title) +
        scale_color_manual(values = c("#1F78B4", "#06A858", "#E31A1C")) +
        coord_flip()
}
