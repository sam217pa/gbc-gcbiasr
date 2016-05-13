#' Conversion tracts of Acinetobacter baylyi
#'
#' A dataset containing the pairwise alignment of \emph{Acinetobacter baylyi}
#' when naturally transformed with a synthetic gene. This synthetic gene
#' introduces SNP that are corrected by MMR. It has been shown that the
#' eucaryotes'MMR is biased towards the conversion to GC. We wanted to look at
#' the conversion bias in bacteria by using this model.
#'
#' @format A data frame with 281,414 rows and 15 variables.
#' \describe{
#'   \item{name}{The experimental name of the transformed clone.}
#'   \item{refb}{The base on the reference sequence}
#'   \item{snpb}{The base on the giver sequence}
#'   \item{expb}{The sequenced base}
#'   \item{refp}{The position on the reference sequence}
#'   \item{cons}{A code giving the sens of the conversion event. \code{.} when the three bases are the same, \code{x} when the base has been converted, \code{X} when the base has been conserved, \code{N} otherwise.}
#'   \item{qual}{The basecall quality}
#'   \item{mutant}{The type of the mutant}
#'   \item{lastmp}{The position of the last marker}
#'   \item{switchp}{The position of the last marker that corresponds to the giver}
#'   \item{switchb}{The base at the breakpoint}
#'   \item{sens}{The sens of the conversion event, W to S (WS) or S to W (SW). }
#'   \item{inconv}{Logical, if the base is in the converted region.}
#'   \item{isrestor}{Logical, if the base is in the converted region, but corresponds to the wild type allele.}
#'   ...
#' }
#' @source \url{http://www.diamondse.info/}
"conversion_tract"
