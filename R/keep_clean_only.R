#' Discard the dirty sequences
#'
#' Discard sequences that are of poor quality, \emph{ie} of a median quality
#' inferior to \code{qual} (40 per default).
#'
#' @title keep_clean_only
#' @param data a dataset, usually the output of \code{read_phruscle}
#' @param qual_ the median quality to discriminate poor quality sequences.
#' @return a data frame filtered off its dirty seq
#' @author Samuel Barreto
#' @import dplyr
#' @importFrom assertthat assert_that
#'
#' @export

keep_clean_only <- function(data, qual_ = 40)
{
    assert_that(
        is.data.frame(data),
        is.numeric(qual_)
    )

    get_dirty_seq <- function(data) {
        dirty_seq <- group_by(data, name) %>%
            summarise(median = median(qual)) %>%
            filter(median < qual_) %>%
            select(name) %>% unlist()
        dirty_marker <- group_by(data, name) %>%
            filter(cons == "x" | cons == "X") %>%
            summarise(median = median(qual)) %>%
            filter(median < qual_) %>%
            select(name) %>% unlist()

        union(dirty_seq, dirty_marker)
    }

    data %>%
        filter(!(name %in% get_dirty_seq(data)))
}
