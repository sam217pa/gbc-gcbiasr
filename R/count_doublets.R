convert_seq_to_char <- function(data) {
    data %>%
        mutate(expb = as.character(expb), snpb = as.character(snpb),
               refb = as.character(refb))
}

#' Make genotype of the sequence
#'
#' Returns a dataframe with the name of the sequence, the type of the mutant,
#' and the total genotype at the markers.
#' @param data a dataframe output of \code{read_phruscle}
#' @param mutant_ a vector of mutant to filter by.
#' @param qual_ the quality bellow which bases are converted to N.
#' @importFrom assertthat assert_that is.string is.count is.flag
#' @import dplyr
#' @export

make_genotype <- function(data
                         ,mutant_ = c("ws", "sw", "s", "w")
                         ,qual_ = 30
                         ,clean=FALSE)
{
    assert_that(
        is.data.frame(data)
       ,is.vector(mutant_) & is.character(mutant_) | is.string(mutant_)
       ,is.count(qual_)
       ,is.flag(clean)
    )

    weak_or_strong <- function(base) {
        is_w <- function(base) { ifelse(base == "A" | base == "T", TRUE, FALSE)}
        is_s <- function(base) { ifelse(base == "C" | base == "G", TRUE, FALSE)}

        if (is_w(base)) "W"
        else if (is_s(base)) "S"
        else if (base == "N") "N"
        else stop(base, " is not a DNA base")
    }

    if (clean) data <- keep_clean_only(data)

    data %>%
        convert_seq_to_char() %>%
        filter(cons == "x" | cons == "X") %>% # ne garde que les positions de SNP
        filter(mutant %in% mutant_) %>%       # filtre selon le mutant
        mutate(isrestor = ifelse(is.na(isrestor) | !isrestor, FALSE, TRUE)) %>%
        rowwise() %>%
        mutate(expb = ifelse(qual < qual_, "N", expb )) %>% # convertit les bases d'une qualité inférieure à qual_ en N
        mutate(genotype = weak_or_strong(expb)) %>% # détermine le génotype de la base séquencée
        ## print()
        mutate(genotype = ifelse(isrestor == TRUE, tolower(genotype), genotype)) %>%
        ungroup() %>% group_by(name, mutant) %>%
        summarise(genotype = toString(genotype) %>% gsub(", ", "", .))  %>% # pretty print
        ungroup()
}

#' La fonction suivante permet de compter le nombre d'anomalies, autrement dit
#' de SNP inattendus dans les traces de conversion. Renvoit une `tbl_df`
#'
#' @param data : le jeu de donné snp
#' @param kmer : le kmer à recenser dans le génotype.

count_kmer <- function(data, kmer) {

    ## closure around find_xxx, which count the number of occurence of the kmer in
    ## the string.
    count_anomaly <- function(string, kmer_) {
        find_xxx <- function(string) {
            function(kmer_) {
                XxX <- gregexpr(kmer_, string) %>% unlist() # compte le nombre d'occurences du kmer
                if   (XxX[1] > 0) length(XxX)
                ## if gregexpr find something, give its length (as integer)
                else 0L # 0 otherwise.
            }
        }

        find_xxx(toupper(string))(kmer_) # count anomaly uses find_xxx to count kmer in string.
    }

    ## applique la fonction count_anomaly par ligne, et ne conserve que les
    ## séquences qui montrent une ou plus d'anomalies.

    data %>%
        rowwise() %>%
        mutate(anom = count_anomaly(genotype, kmer)) %>%
        filter(anom > 0)
}

#' Count last SNP pattern
#'
#' Per type of mutant, count the last SNP of the gene conversion tract,
#' and determine its type, either W (AT) or S (GC).
#'
#' @param data a dataframe results of \code{read_phruscle}
#' @return a table
#' @author Samuel Barreto
#' @import dplyr
#' @import ggplot2
#' @importFrom assertthat is.flag
#' @export

count_last_snp <- function(data)
{
    if (!(is.data.frame(data)))
        stop(data, " is not a dataframe")

    data <-
        data %>%
        filter(switchp != lastmp, refp == switchp, cons == "x" | cons == "X") %>%
        ## print()
        ## filter(is.na(sens))
        group_by(mutant, sens) %>%
        summarise(count = n()) %>%
        ungroup() %>%
        mutate(mutant = factor(mutant, labels = c("S", "W", "SW", "WS"))
              ,sens = factor(sens, levels = c("ws", "sw"), labels = c("CG", "AT")))

    class(data) <- c("lastsnp", class(data))

    data
}


#' @title Count restoration of wild type haplotype
#' @description
#' This function count the number of event of wild type haplotype restoration.
#'
#' @param data a data frame of snp calling.
#' @param quality quality bellow which bases are not to be trusted
#' @importFrom assertthat is.flag is.count
#' @import dplyr
#'
#' @export
count_restor <- function(data, quality = 40)
{
    if (!(is.data.frame(data)))
        stop(data, " must be a dataframe")
    if (!is.count(quality))
        stop(quality, " is not a valid quality")

    data <-
        data %>%
        filter(isrestor, qual >= quality) %>%
        ## mutate(mutant = factor(mutant, labels = c("S", "W", "SW", "WS"))) %>%
        group_by(mutant, sens) %>%
        summarise(count = n()) %>%
        ungroup() %>%
        mutate(sens = factor(sens, labels = c("CG", "AT")))

    class(data) <- c("lastsnp", class(data))
    data
}

#' @export
print.lastsnp <- function(x, ...)
{
    cat(knitr::kable(x, align = "c"), sep = "\n")
}

#' @export
plot.lastsnp <- function(x, ...)
{
    ggplot(x, aes(x = mutant, y = count, color = sens )) +
        geom_point() +
        coord_flip() +
        scale_color_solarized(
            labels = c("AT", "CG"),
            guide = guide_legend(title = "SNP au point\nde bascule")) +
        scale_y_continuous(breaks = extended_range_breaks()(x$count)) +
        labs(x = "Donneur", y = "") +
        theme(legend.position = "right")
}
