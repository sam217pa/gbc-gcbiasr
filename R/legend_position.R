#' Easily place the legend on the graph.
#'
#' If both x and y are null, place it at the bottom by default.
#' @title legend_position
#' @param x the x coordinate. must be between 0 and 1
#' @param y idem.
#' @return a theme element.
#' @author Samuel Barreto
#' @importFrom ggplot2 theme
#' @importFrom scales alpha
#' @export

legend_position <- function(x=NULL, y=NULL) {
    custom_legend <- theme(
        legend.margin = unit(-0.3,"lines"),
        legend.key = element_rect(fill = scales::alpha("gray", 0),
                                  colour = scales::alpha("gray", 0)))
    if (is.null(x) & is.null(y)) {
        custom_legend + theme(legend.position = "bottom")
    } else {
        custom_legend + theme(legend.position = c(x, y))
    }
}
