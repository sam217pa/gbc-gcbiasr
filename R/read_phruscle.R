#' Read the csv dataset in ``path'' and return it as a clean dataset.
#'
#' @title read_phruscle
#' @param path the path to the csv file.
#' @return a tbl_df
#' @author Samuel Barreto
#' @importFrom assertthat assert_that has_extension is.string
#' @importFrom readr read_csv
#' @import dplyr
#'
#' @export

read_phruscle <- function(path=NULL, mutant_levels = NULL)
{
    ## test
    assert_that(
        file.exists(path),
        has_extension(path = path, ext = "csv"),
        is.null(mutant_levels) | is.vector(mutant_levels) | is.string(mutant_levels)
    )

    ## default parameters
    if (is.null(mutant_levels)) mutant_levels <- c("s", "w", "sw", "ws")

    ## read data, set default column types
    snp <- read_csv(path, col_types = "ccciciccid", trim_ws = TRUE)

    ## set the correct base levels
    set_base_level <- function(x) factor(x, levels = c("G", "C", "T", "A", "N"))

    snp <- mutate(
        snp,
        name = gsub("-1073.+$", "", name),
        refb = set_base_level(refb),
        snpb = set_base_level(snpb),
        seqb = set_base_level(seqb)
    ) %>%
        select(name = name, refb, snpb, expb = seqb, refp, expp = seqp, cons,
               qual, base)

    if (FALSE %in% (toupper(snp$base) == snp$expb)) {
        warning("Bases are not the same in the phred file and the parsed alignment.",
                "You should check them, they're usually the sign of",
                "alignment artifacts")
    } else {
        snp <- select(snp, -base)
    }


    ## add a column which set the mutant type
    find_mutant <- function(data)
    {
        get_mutant_type <- function(name)
        {
            if      (grepl("ws", name)) "ws"
            else if (grepl("sw", name)) "sw"
            else if (grepl("W", name )) "w"
            else if (grepl("w", name )) "w"
            else if (grepl("S", name )) "s"
            else if (grepl("s", name )) "s"
            else ""
        }

        ## neat little trick to reduce time of rowwise application of find_mutant.
        ## get the mutant's type, either ws, sw, w or s.
        data %>%
            group_by(name) %>%
            summarise(id = unique(name)) %>%
            rowwise() %>%
            mutate(mutant = get_mutant_type(id)) %>%
            mutate(mutant = factor(mutant, levels = mutant_levels)) %>%
            inner_join(data, .) %>%
            select(-id) %>%
            filter(mutant != "")
    }

    ## add a column which sets the lastmp, last marker position, starting from
    ## the conversion side.
    find_lastmp <- function(data)
    {
        data %>% group_by(name) %>%
            filter(cons == "x" | cons == "X") %>%
            summarise(lastmp = min(refp)) %>%
            inner_join(data, .)
    }

    ## add a column which sets the the switchp, switch position.
    find_switchp <- function(data)
    {
        max_refp <- max(data$refp)
        switch_pos_table <- rbind(
            ## deals with cases where there is a conversion tract.
            data %>%
            group_by(name) %>%
            filter(cons == "x") %>%
            summarise(switchp = min(refp)),
            ## deals with cases where there is no conversion tract.
            data %>%
            group_by(name) %>%
            filter(!("x" %in% cons)) %>%
            summarise(switchp = max_refp)
        )
        inner_join(data, switch_pos_table)
    }

    ## add a column which sets the base at the switch position per sequence
    find_switchb <- function(data) {
        data %>%
            filter(cons == "x", refp == switchp) %>%
            mutate(switchb = expb) %>%
            select(name, switchb) %>%
            left_join(data, .)
    }

    ## détermine si un snp est dans le sens W to S ou S to W
    find_polarity <- function(data)
    {

        strong_or_weak <- function(base)
        {
            is_w <- function(base) { ifelse(base == "A" | base == "T", TRUE, FALSE)}
            is_s <- function(base) { ifelse(base == "C" | base == "G", TRUE, FALSE)}

            if      (is_w(base)) "w"
            else if (is_s(base)) "s"
            else
                stop(base, "is not a dna base")
        }

        data %>%
            filter(cons == "x" | cons == "X") %>%
            rowwise() %>%
            mutate(sens = paste0(strong_or_weak(refb),
                                 strong_or_weak(expb))) %>%
            ungroup() %>%
            left_join(data, .)
    }

    find_isinconv <- function(data)
    {
        data %>%
            ## si la position est avant la position de switch.
            mutate(inconv = ifelse(refp >= switchp, TRUE, FALSE)) %>%
            inner_join(data, .)
    }

    find_isrestor <- function(data)
    {
        data %>%
            filter(inconv, cons == "x" | cons == "X", expb == refb) %>%
            mutate(isrestor = TRUE) %>%
            left_join(data, .)
    }

    snp %>%
        find_mutant() %>%
        find_lastmp() %>%
        find_switchp() %>%
        find_switchb() %>%
        find_polarity() %>%
        find_isinconv() %>%
        find_isrestor()
}
