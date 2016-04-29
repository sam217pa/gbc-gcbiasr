##' Plot the alignment of the markers positions
##'
#' @title plot_align
#' @param data a data.frame. it must be the output of \code{read_phruscle}.
#' @param mutant_ the mutant to filter by.
#' @param plot_title chose a title. if NULL, guess from the mutant parameter
#' @param quality min quality to consider a restoration as such
#' @return a ggplot output of the alignment
#' @author Samuel Barreto
#' @import dplyr
#' @import ggplot2
#' @import cowplot
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggthemes extended_range_breaks
#' @importFrom assertthat assert_that is.string is.flag
#' @export
plot_align <- function(data, mutant_, plot_title=NULL, quality = 30, color=TRUE)
{
    assert_that(
        is.data.frame(data)
       ,is.string(mutant_)
       ,is.null(plot_title) | is.string(plot_title)
       ,is.flag(color)
    )

    if (is.null(plot_title))
        plot_title <- paste("Alignement pour la manip", toupper(mutant_))

    ## default color
    rec.col  <- "#E45052" # red paler than Set1
    don.col  <- "#72BCE4"
    rest.col <- brewer.pal(n = 4, "Spectral")[3]
    ## gc.col  <- brewer.pal(n = 4, "Spectral")[4]

    sort_by_tract_length <- function(data, mut)
    {
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
        filter(cons == "x" | cons == "X") %>%
        rowwise() %>%
        mutate(expb = ifelse(expb %in% c("A", "T"), "w", "s"))

    first_snp <- filter(data, cons == "X") %>% summarise(min = min(refp))
    last_snp <- filter(data, cons == "x") %>% summarise(max = max(refp))
    num_of_seq <- n_distinct(data$name)

    ## default output theme.
    set_gcbiasr_theme()

    align_plot <- ggplot(data = data, aes(x = refp, y = name)) +
        geom_point(aes(color = inconv, size = qual, alpha=qual)) +
        ## représente la séquence du donneur
        geom_text(aes(label = snpb, x = refp, y = -7), color = don.col,
                  vjust = -2, size = 4, family = "Ubuntu Light") +
        geom_text(aes(label = refb, x = refp, y = num_of_seq + 10)
                 ,color = rec.col, vjust = 4, size = 4
                 ,family = "Ubuntu Light") +
        scale_alpha( range=c(1/5, 0.8), guide=FALSE ) +
        scale_size(range = c(0.5, 2), breaks = c(10, 50),
                   labels = c("Faible", "Forte")) +
        scale_x_continuous(breaks = extended_range_breaks()(data$refp)) +
        ## du premier au dernier snp
        coord_cartesian(xlim = c(first_snp, last_snp)) +
        labs(x = "", y = "", size = "Qualité" ,title = plot_title
            ,color = "Haplotype" ,shape = "Restauration\nvers…") +
        guides(
            ## change le symbol et la taille de la légende pour les couleurs
            colour = guide_legend(override.aes = list(shape = 20, size = 3))
            ## change la couleur par défault de la légende de taille
           ,size  = guide_legend(override.aes = list(color = don.col))) +
        theme(legend.margin = unit(0,"lines"),
              panel.grid.major.y = element_line(size = 0.1, linetype = "dotted"),
              legend.position = "right",
              legend.justification = c(0, 1),
              legend.text = element_text(size = 6))

    ## représente les cas de traces de conversion complexes
    align_plot <- align_plot +
        geom_point(data = filter(data, isrestor == TRUE, qual > quality)
                  ,aes(x = refp + 12 ,shape = expb)
                  ,size = 3) +
        scale_shape_manual(breaks = c("w", "s")
                          ,values = c("S", "W")
                          ,labels = c("AT", "GC"))
    ## }

    if (color) {
        align_plot <- align_plot +
            ## représente les cas complexes
            geom_point(data = filter(data, isrestor == TRUE, qual > quality),
                       aes(size = qual),
                       color = rec.col) +
            scale_color_manual(
                limits = c("TRUE", "FALSE")
               ,values = c(don.col, rec.col)
               ,labels = c("Donneur", "Receveur"))

    } else {
        align_plot <- align_plot +
            ## représente les cas complexes
            geom_point(data = filter(data, isrestor == TRUE, qual > quality)
                      ,aes(size = qual)
                      ,color = rec.col) +
            scale_color_manual(
                limits = c("TRUE", "FALSE")
               ,values = c("black", "gray")
               ,labels = c("Donneur", "Receveur"))

    }

    ggdraw(align_plot) +
        draw_label(label = "Donneur", x = 0.07, y = 0.1,
                   fontfamily = "Ubuntu Light", size = 8, colour = don.col) +
        draw_label(label = "Receveur", x = 0.07, y = 0.89, colour = rec.col,
                   fontfamily = "Ubuntu Light", size = 8)
    ## align_plot
}
