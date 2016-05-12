##' Plot the alignment of the markers positions
##'
#' @title plot_align
#' @param data a data.frame. it must be the output of \code{read_phruscle}.
#' @param mutant_ the mutant to filter by.
#' @param plot_title chose a title. if NULL, guess from the mutant parameter
#' @param quality min quality to consider a restoration as such
#' @param color logical : use color in the conversion tract or black and white.
#'     color by default
#' @param inverse logical : must inverse the x scale or not. inverse by default
#' @param with_lab logical : keep y labels or not ? usually you don't need them.
#' @return a ggplot output of the alignment
#' @author Samuel Barreto
#' @import dplyr
#' @import ggplot2
#' @importFrom cowplot ggdraw draw_label
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggthemes extended_range_breaks
#' @importFrom assertthat assert_that is.string is.flag
#' @export
plot_align <- function(data, mutant_, plot_title=NULL, quality = 30
                      ,color=TRUE, inverse=TRUE, with_lab=FALSE)
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
        mutate(expb = ifelse(expb %in% c("A", "T"), "w", "s")) %>%
        ungroup()

    first_snp <- filter(data, cons == "X") %>% summarise(min = min(refp))
    last_snp <- filter(data, cons == "x") %>% summarise(max = max(refp))
    num_of_seq <- n_distinct(data$name)
    label_table <- summarise(
        group_by(data, refp)
       ,don_label = unique(snpb)
       ,rec_label = unique(refb)
       ,ref_pos   = unique(refp)
    )
    ## default output theme.
    set_gcbiasr_theme()

    align_plot <- ggplot(data = data, aes(x = refp, y = name)) +
        geom_point(aes(color = inconv, size = qual, alpha=qual)) +
        ## représente la séquence du donneur
        geom_text(data = label_table, aes(label = don_label, x = ref_pos, y = -1), color = don.col,
                  size = 3, family = "Ubuntu") +
        geom_text(data = label_table, aes(label = rec_label, x = ref_pos, y = num_of_seq + 2)
                 ,color = rec.col, size = 3 ,family = "Ubuntu") +
        scale_alpha( range=c(1/5, 0.8), guide=FALSE ) +
        scale_size(range = c(0.5, 2), breaks = c(10, 50),
                   labels = c("Faible", "Forte")) +
        ## scale_y_discrete()
        ## du premier au dernier snp
        coord_cartesian(ylim = c(-3, num_of_seq + 3)) +
        labs(x = "", y = "", size = "Qualité" ,title = plot_title
            ,color = "Haplotype" ,shape = "Restauration\nvers…") +
        guides(
            ## change le symbol et la taille de la légende pour les couleurs
            colour = guide_legend(override.aes = list(shape = 20, size = 3))
            ## change la couleur par défault de la légende de taille
           ,size  = guide_legend(override.aes = list(color = don.col))) +
        theme(plot.title = element_text(size = 10),
              legend.margin = unit(0,"lines"),
              panel.grid.major.y = element_line(size = 0.1, linetype = "dotted"),
              legend.position = "right",
              legend.justification = c(0, 1),
              legend.text = element_text(size = 6))

    ## représente les cas de traces de conversion complexes {
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

    ## if inverse is TRUE, reverse the x scale. place the origin of conversion to 0.
    ## FIXME pour une raison que j'ignore, ça ne marche pas.
    if (inverse) {
        align_plot <- align_plot +
            scale_x_continuous(breaks = extended_range_breaks()(data$refp)
                              ,trans = "reverse")
    } else {
        align_plot <- align_plot +
            scale_x_continuous(breaks = extended_range_breaks()(data$refp))
    }

    if (!with_lab) {
        align_plot <- align_plot +
            theme(axis.text.y = element_blank())
    }

    ## Je me suis dit que pour les graphiques finaux on n'a pas besoin forcément des
    ## labels. On peut donc dans le theme supprimer les étiquettes de l'axe des y.

    ## ggdraw(align_plot) +
    ##     draw_label(label = "Donneur", x = 0.07, y = 0.1,
    ##                fontfamily = "Ubuntu Light", size = 8, colour = don.col) +
    ##     draw_label(label = "Receveur", x = 0.07, y = 0.89, colour = rec.col,
    ##                fontfamily = "Ubuntu Light", size = 8)
    align_plot
}
