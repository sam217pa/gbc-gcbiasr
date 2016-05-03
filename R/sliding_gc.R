
#' Measure GC content by sliding window
#'
#' Measure the GC content of a sequence by a sliding window of a given size.
#'
#' @details
#' It uses the GC function of the seqinr package, and the rollapply from the
#' zoo package.
#'
#' @param sequence a character vector of the sequence
#' @importFrom zoo zoo rollapply
#' @import seqinr 
#' @import dplyr
#'
#' @export

sliding_gc <- function(sequence, window_width = 100, slide_by = 1)
{
    if (!(file.exists(sequence)))
        stop(sequence, " does not exists at the current adress")
    if (!(is.count(window_width)))
        stop("Window width must be a positive integer")
    if (!(is.count(slide_by)))
        stop("slide_by must be a positive integer")

    fasta <- seqinr::read.fasta(sequence)

    seq.zoo <- zoo::zoo(seqinr::getSequence(fasta)[[1]])

    seq.gc.sliding <- zoo::rollapply(
	    seq.zoo
	    ,width = window_width
	    ,by = slide_by
	    ,FUN = GC
	    ,align = "center" 
	    ,fill = NA
	    )

    sequence.gc.data <-
	    data_frame(seq.gc.sliding) %>%
	    add_rownames() %>%
	    select(id = rowname, GC = seq.gc.sliding) %>%
	    mutate(id = as.numeric(id))

    class(sequence.gc.data) <- c("gcslide", class(sequence.gc.data))

    sequence.gc.data
}


#' @import ggplot2
#'
#' @export

plot.gcslide <- function(x, plot_title=NULL, wwidth=NULL, ...)
{
    wwidth <- attr(x, "window_width")

    if (is.null(plot_title))
        ## plot_title <- paste0("GC moyen\n", "Fenêtre : ", wwidth, "pb")
        plot_title <- "GC moyen"

    ggplot(data = x, aes(x = id, y = GC)) +
        geom_line(color = "gray") +
        ## geom_point(color = "black", size = 1, alpha = 0.2) +
        ## geom_smooth(span = span, se = FALSE) +
        labs(x = "Position\nsur la référénce", y = "", title = plot_title)
}
