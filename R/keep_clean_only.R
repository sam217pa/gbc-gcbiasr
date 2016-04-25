##' Discard the dirty sequence
##'
##' @title keep_clean_only
##' @param data a dataset, usually the output of \code{read_phruscle}
##' @param qual_ the median quality to discriminate poor quality sequences.
##' @return a data frame filtered off its dirty seq
##' @author Samuel Barreto
##' @import dplyr
##' @export
keep_clean_only <- function(data, qual_ = 40) {

    get_dirty_seq <- function(data)
        group_by(data, name) %>%
            summarise(median = median(qual)) %>%
            filter(median < qual_) %>%
            select(name) %>% unlist()

    data %>%
        filter(!(name %in% get_dirty_seq(data)))
}
