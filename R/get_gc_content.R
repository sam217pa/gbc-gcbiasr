
#' compute the gc content per sequence
#'
#' @param data a dataframe output of \code{read_phruscle}
#' @param mut the mutant to filter by.
#' @importFrom seqinr GC
#' @importFrom assertthat assert_that is.string
#' @importFrom tidyr gather
#' @import dplyr
#'
#' @export

get_gc_content <- function(data, mut) {

    assert_that(
        is.data.frame(data),
        is.null(mut) | is.vector(mut) | is.string(mut)
    )

    if (is.null(mut)) mut <- c("ws", "sw", "w", "s")

    gc_content <- data %>%
        ## omit the rare cases of N in sequences
        filter(data, !(is.na(refb)), !(is.na(snpb)), !(is.na(expb))) %>%
        mutate(snpb = as.character(snpb),
               expb = as.character(expb),
               refb = as.character(refb)) %>%
        filter(mutant %in% mut) %>%
        group_by(name, mutant) %>%
        ## compute gc content of the respective sequence
        summarise(exp = seqinr::GC(expb),
                  snp = seqinr::GC(snpb),
                  ref = seqinr::GC(refb))

    ## convert the data set to long form.
    gc_content %>%
        gather("type", "GC", 3:5) %>%
        mutate(type = factor(
                   type, levels = c("snp", "exp", "ref"),
                   labels = c("Donneur", "Recombinant", "Receveur")))
}
