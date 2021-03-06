##' Set the default theme for ggplot2 plots.
##'
##' Inspired both by the fivethirtyeight theme and the solarized theme of
##' ggthemes.
##'
##' @title set_gcbiasr_theme
##' @return the default theme for ggplot2 plots
##' @author Samuel Barreto
##' @import ggplot2
##' @importFrom RColorBrewer brewer.pal
##' @export
set_gcbiasr_theme <- function() {

    fte_theme <- function() {
        ## inspiré de http://minimaxir.com/2015/02/ggplot-tutorial/
        ## Generate the colors for the chart procedurally with RColorBrewer
        palette <- brewer.pal("Greys", n=9)
        plot.background = "white"
        color.background = "#F8F6EE"
        color.strip.background = "white"
        ## alpha(color.background, 1)
        color.grid.major = palette[3]
        color.axis.text = palette[6]
        color.axis.title = palette[7]
        color.title = palette[9]
        ## Begin construction of chart
        theme_bw(base_size=9) +
            ## Set the entire chart region to a light gray color
            theme(
                text = element_text(family = "Ubuntu"),
                panel.background=element_rect(
                    fill=alpha(color.background, 1), color=color.background),
                plot.background=element_rect(
                    fill=plot.background, color = plot.background),
                plot.background = element_blank(),
                panel.border=element_rect(color=color.background),
                panel.grid.major=element_line(
                    color=color.grid.major,size=.25, linetype = "dotted"),
                panel.grid.minor=element_blank(),
                axis.ticks=element_blank(),
                strip.text.y = element_text(angle = 0, size = 8),
                strip.background = element_rect(fill = color.strip.background,
                                                color = color.strip.background),
                legend.position="none", # Format the legend, but hide by default
                legend.background = element_blank(),
                legend.text = element_text(size=6,color=color.axis.title),
                legend.title = element_text(size=8,color=color.axis.title),
                legend.key = element_rect(
                    fill = scales::alpha("gray", 0),
                    colour = scales::alpha("gray", 0)),
                ## Set title and axis labels, and format these and tick marks
                plot.title = element_text(
                    color = color.title, size = 14, face = "bold", vjust = 1.25,
                    hjust = 0, family = "Ubuntu"),
                axis.text.x = element_text(size = 7,color = color.axis.text),
                axis.text.y = element_text(size = 7,color = color.axis.text),
                axis.title.x = element_text(
                    size = 8, color = color.axis.title, vjust = 0, hjust = 0.9),
                axis.title.y = element_text(
                    size = 8, color = color.axis.title, vjust = 0.9, angle = 0 ),
                ## Plot margins
                plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
            )
    }

    theme_set(theme_bw() + fte_theme())
}
