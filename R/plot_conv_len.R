#' Plot the distribution of the conversion tracts
#' 
#' This function takes a dataset as output of \code{read_phruscle} and 
#' plot the distribution of the conversion tracts.
#' @param data the output of \code{read_phruscle}
#' @return a ggplot2 histogram of the distribution of conversion tracts.
#' @import dplyr
#' @import ggplot2
#' @export

plot_conv_len <- function(data)
{
        weak_or_strong <- function(base) {

                is_w_s <- function(base)
                        function(motif) ifelse(base %in% motif, T, F)

                if      (is_w_s(base)(c("A", "T"))) "w"
                else if (is_w_s(base)(c("C", "G"))) "s"
                else "N"
        }


        data %>%
                group_by(mutant, name ) %>%
                summarise(switchp = unique(switchp), switchb = unique(switchb)) %>%
                group_by(mutant, switchp, switchb) %>%
                summarise(distri = n()) %>%
                rowwise() %>%
                mutate(switchb = weak_or_strong(switchb)) %>%
                ggplot(aes(x = switchp, y = distri, fill = switchb)) +
                geom_bar(stat = "identity") +
                geom_point(aes(color = switchb)) +
                #stat_summary(fun.y = mean, geom = "point", color = "black", shape = "¦",
                             #size = 6) +
                facet_grid(mutant ~. ) +
                scale_fill_brewer(palette = "Set1", guide = FALSE) +
                scale_color_brewer(palette = "Set1", labels = c("N", "GC", "AT")) +
                labs(x = "Position de changement d'haplotype\nentre donneur et receveur",
                     y = "",
                     color = "Base au\npoint de bascule") +
                legend_position()
}
