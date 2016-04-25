##' Plot the alignment of the markers positions
##'
##' @title plot_align
##' @param data a data.frame. it must be the output of \code{read_phruscle}.
##' @param mutant_ the mutant to filter by.
##' @return a ggplot output of the alignment
##' @author Samuel Barreto
##' @import dplyr
##' @import ggplot2
##' @importFrom RColorBrewer brewer.pal
##' @importFrom ggthemes extended_range_breaks
##' @export
plot_align <- function(data, mutant_) {
    ## default color
    red    <- brewer.pal(n = 4, "Set1")[1]
    blue   <- brewer.pal(n = 4, "Set1")[2]
    green  <- brewer.pal(n = 4, "Set1")[3]
    violet <- brewer.pal(n = 4, "Set1")[4]

    sort_by_tract_length <- function(data, mut) {
        ## combine les données pour les clones avec conversion et les clones
        ## sans, pour lesquels on ne peut pas déterminer de longueur de trace de
        ## conversion.
        sorted_seq <- data %>%
            group_by(name) %>%
            summarise(len = unique(switchp)) %>%
            ## trie par longueur de trace
            arrange(len) %>%
            ## convertit en vecteur
            select(name) %>% unlist() %>% as.vector()

        data %>%
            mutate(name = factor(name, levels = sorted_seq))
    }

    data <- data %>%
        filter(mutant %in% mutant_) %>%
        keep_clean_only() %>%
        sort_by_tract_length() %>%
        filter(cons == "x" | cons == "X")

    first_snp <- filter(data, cons == "X") %>% summarise(min = min(refp))
    last_snp <- filter(data, cons == "x") %>% summarise(max = max(refp))
    num_of_seq <- n_distinct(data$name)

    ## default output theme.
    set_gcbiasr_theme()

    ggplot(data = data, aes(x = refp, y = name)) +
        geom_point(aes(color = inconv, size = qual, alpha=qual)) +
        ## représente les cas complexes
        geom_point(data = filter(data, isrestor == TRUE),
                   ## aes(alpha = qual), color = green, size = 5) +
                   aes(alpha = qual), color = red, size = 3) +
        ## représente la séquence du donneur
        geom_text(aes(label = snpb, x = refp, y = -7),
                  color = "black",
                  vjust = -2,
                  size = 4,
                  family = "Ubuntu Light") +
        ## annotate("text", x = as.numeric(first_snp) - 20, y = -7, vjust = -2,
        ##          label = "Donneur : ", color = blue, size = 4 ) +
        ## ## représente la séquence du receveur.
        geom_text(aes(label = refb, x = refp,
                      y = num_of_seq + 10),
                  color = "grey",
                  vjust = 4,
                  size = 4,
                  family = "Ubuntu Light") +
        ## annotate("text", x = as.numeric(first_snp) - 20,
        ##          y = num_of_seq + 10, vjust = 4,
        ##          label = "Receveur : ", color = red, size = 4 ) +
        ## j'ai tenté de représenter les séquences par des segments de couleur,
        ## mais je trouve que ça perturbe la lecture. On a un effet moiré pas du
        ## tout attendu, qui n'apporte rien.
        ##
        ## geom_segment(data =
        ##                  filter(data, cons == "x") %>%
        ##                  group_by(name) %>%
        ##                  summarise(max = max(refp), min = min(refp)),
        ##              aes(x = max, xend = min, y = name, yend = name),
        ##              color = blue, alpha = 0.2) +
        ## geom_segment(data =
        ##                  filter(data, cons == "X") %>%
        ##                  group_by(name) %>%
        ##                  summarise(max = max(refp), min = min(refp)),
        ##              aes(x = max, xend = min, y = name, yend = name),
        ##              color = red, alpha = 0.2) +
        scale_color_manual(values = c("gray", "black"), guide = FALSE) +
        ## scale_color_brewer(palette = "Set1", guide = "none") +
        scale_alpha( range=c(1/5, 0.8), guide=FALSE ) +
        scale_size(range = c(0.5, 2), breaks = c(10, 50),
                   labels = c("Faible", "Forte")) +
        scale_x_continuous(breaks = extended_range_breaks()(data$refp)) +
        coord_cartesian(xlim = c(first_snp, last_snp)) +
        labs(x = "", y = "", size = "Qualité",
             title = paste("Alignement pour la manip", toupper(mutant_))) +
        theme(legend.direction = "horizontal",
              legend.margin = unit(0,"lines"),
              panel.grid.major.y = element_line(size = 0.1, linetype = "dotted")) +
        legend_position(0.8, 0.98)
}

##' Plot the quality of the sequence
##'
##' @title plot_qual
##' @param data the output of \code{run_phruscle}.
##' @param mutant_ the mutant to filter by.
##' @return a ggplot2 plot.
##' @author Samuel Barreto
##' @importFrom RColorBrewer brewer.pal
##' @import dplyr
##' @import ggplot2
##' @export

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


##' Plot the alignment and the quality of the sequence juxtaposed
##'
##'
##' @title plot_align_qual
##' @param data the output of \code{run_phruscle}.
##' @param mutant_ the mutant to filter by.
##' @return a ggplot2 plot.
##' @author Samuel Barreto
##' @importFrom cowplot plot_grid
##' @export

plot_align_qual <- function(data, mutant_) {
    plot_grid(
        plot_align(data, mutant_),
        plot_qual(data, mutant_),
        ncol = 1, rel_heights = c(8.5, 1.5),
        align = 'v', labels = c("  1", "  2")
    )
}
