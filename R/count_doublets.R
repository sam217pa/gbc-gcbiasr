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

make_genotype <- function(data, mutant_ = c("ws", "sw", "s", "w"), qual_ = 30, clean=FALSE) {

    assert_that(
        is.data.frame(data),
        is.vector(mutant_) & is.character(mutant_) | is.string(mutant_),
        is.count(qual_),
        is.flag(clean)
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

#' Count pattern in genotype
#'
#' Count the occurence of \code{kmer} in the genotype of mut.
#' @param data a dataframe results of \code{read_phruscle}
#' @param mut the mutant to filter by
#' @param kmer the kmer to count.
#' @return a table
#' @author Samuel Barreto
#' @import dplyr
#' @importFrom assertthat is.string is.flag
#'
#' @export

count_motif <- function(data, mut, kmer, clean=TRUE)
{
    valid_kmer <- c("WWW", "SWW", "WSW", "SSW", "WWS", "SWS", "WSS", "SSS",
                    "WW", "SW", "WS", "SS")

    if (!is.data.frame(data))
        stop(data, "must be a dataframe.")
    if (!(is.vector(mut) & is.character(mut) | is.string(mut)))
        stop(mut, " must be a vector of character or a single string.")
    if (!is.string(kmer))
        stop(kmer, " must be a countable motif.")
    if (!(kmer %in% valid_kmer))
        stop(kmer, " is not a valid kmer to search for.")
    if (!is.flag(clean))
        stop(clean, " should be logical.")

    num_of_seq <- n_distinct(filter(data, mutant %in% mut)$name)

    ## compte les kmer.
    kmer_counter <- function(kmer)
    {
        data %>% make_genotype(mut, clean = clean) %>% count_kmer(kmer) %>%
            ungroup() %>%
            { if (nrow( . ) > 0) select( . , anom) %>% sum() else 0}
    }

    if (nchar(kmer) == 3) {
        kmer_counter(kmer)
    } else {
        ## do something else. must take into account the restoration. could
        ## change the count of doublets
        kmer <- unlist(strsplit(kmer, ""))
        kmer_xy <- gsub(", ", "", toString(c(kmer[1], kmer[2])))
        kmer_xyy <- gsub(", ", "", toString(c(kmer[1], kmer[2], kmer[2])))
        kmer_yxx <- gsub(", ", "", toString(c(kmer[2], kmer[1], kmer[1])))

        if (kmer[1] == kmer[2]) {
            total_count <- kmer_counter(kmer_xy) -
                kmer_counter(kmer_xyy)
        } else {
            total_count <- kmer_counter(kmer_xy) -
                kmer_counter(kmer_xyy) - kmer_counter(kmer_yxx)
        }

        if (!(total_count < num_of_seq)) {
            stop(total_count, " is greater than the number of sequences.\n",
                 "Check that kmer and mut are coherent values to look for.")
        } else if (!(total_count > 0)) {
            stop(total_count, " is less than 0\n",
                 "Check that kmer and mut are coherent values to look for.")
        } else {
            total_count
        }
    }
}
