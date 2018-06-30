#' Binding sites for selected restriction enzymes.
#'
#' @format A 695 x 2 data frame with the following columns:
#' \describe{
#'   \item{enzyme}{Restriction enzyme name.}
#'   \item{sites}{Enzyme binding site sequences.
#'       See \code{\link{nucleobase_legend}} for what bases other than `T`, `C`, `A`,
#'       and `G` mean.
#'       Each `/` indicates a cleavage point.
#'       According to NEB (see link below under Source):
#'       "Numbers in parentheses indicate point of cleaveage for non-palindromic
#'       enzymes."
#'       These types of enzymes are not implemented.}
#' }
#' @source \url{https://www.neb.com/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities}
"binding_sites"



#' Legend for the single-letter code of nucleobases indicating restriction sequences.
#'
#' @format A data frame of 28 rows and two columns:
#' \describe{
#'   \item{code}{The letter indicating more than one possible nucleotides.}
#'   \item{nucleotides}{One of the multiple nucleotides the code can refer to.}
#' }
#' @source \url{https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites:_A#Whole_list_navigation}
"nucleobase_legend"




#' Table of evolutionary rates.
#'
#' From Table 1 in Sung et al. (2016).
#'
#'
#' @format A data frame with 15 rows and 8 variables:
#' \describe{
#'   \item{domain}{Either \code{Bacteria} or \code{Eukarya} for what type of organism
#'       the species is.}
#'   \item{species}{Species name.}
#'   \item{Ge}{Effective genome size using only coding DNA.}
#'   \item{Gc_Gnc}{Effective genome size using coding DNA and non-coding DNA that is
#'       under purifying selection.}
#'   \item{indels}{Rate of insertions and deletions (\eqn{\times 10^{-10}}
#'       events per site per generation).}
#'   \item{subs}{Base-substitution mutation rate (\eqn{\times 10^{-10}}
#'       events per site per generation).}
#'   \item{Ne}{Effective population size (\eqn{\times 10^6}).}
#'   \item{theta_s}{Population mutation rate estimated using \eqn{{\theta}_s}.}
#'   \item{pi_s}{Population mutation rate estimated using \eqn{{\pi}_s}.}
#' }
#'
#' @source \url{http://dx.doi.org/10.1534/g3.116.030890}
#'
#' @references
#' Sung, W., M. S. Ackerman, M. M. Dillon, T. G. Platt, C. Fuqua, V. S. Cooper, and
#' M. Lynch. 2016. Evolution of the insertion-deletion mutation rate across the
#' tree of life. \emph{G3: Genes | Genomes | Genetics} \strong{6}:2583–2591.
#'
#'
"evo_rates"


