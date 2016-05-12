#' Plot the distribution of the conversion tracts
#'
#' This function takes a dataset as output of \code{read_phruscle} and
#' plot the distribution of the conversion tracts.
#' @param data the output of \code{read_phruscle}
#' @return a ggplot2 histogram of the distribution of conversion tracts.
#' @import dplyr
#' @import ggplot2
#' @export

breakpoints_distribution <- function(data)
{
    weak_or_strong <- function(base) {

        is_w_s <- function(base)
            function(motif) ifelse(base %in% motif, T, F)

        if      (is_w_s(base)(c("A", "T"))) "w"
        else if (is_w_s(base)(c("C", "G"))) "s"
        else "N"
    }

    # TODO commenter code. rien n'est très clair à la relecture

    distri_snp <- data %>%
        group_by(mutant, name ) %>%
        summarise(switchp = unique(switchp), switchb = unique(switchb)) %>%
        group_by(mutant, switchp, switchb) %>%
        summarise(distri = n()) %>%
        rowwise() %>%
        mutate(switchb = weak_or_strong(switchb)) %>%
        mutate(switchb = factor(switchb, levels = c("s", "w", "N")))

    class(distri_snp) <- c("breakpoints", class(distri_snp))

    distri_snp
}


#' @import ggplot2
#' @export

plot.breakpoints <- function(x, facet=TRUE, title = "", ...)
{
    p <-
        ggplot(data = x, aes(x = switchp, y = distri, fill = switchb)) +
        scale_fill_solarized(guide = FALSE) +
        scale_color_solarized(labels = c("GC", "AT", "N")) +
        labs(x = "Position de changement d'haplotype\nentre donneur et receveur"
            ,y = ""
            ,color = "Base au\npoint de bascule"
            ,title = title) +
        theme(legend.position = "top")

    if (facet) p + facet_grid(mutant ~. ) +
                   geom_bar(stat = "identity") +
                   geom_point(aes(color = switchb))
    else p + geom_bar(stat = "identity"
                     ,position = "stack"
                     ,width = 1)

}
