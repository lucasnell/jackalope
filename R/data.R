#' Binding sites for selected restriction enzymes.
#'
#' @format A list of length 11. For each element in the list...
#' \describe{
#'   \item{name}{The element's name is the restriction enzyme's name.
#'       The enzymes present in this list are
#'       \emph{AclI}, \emph{ApeKI}, \emph{AscI}, \emph{BspEI}, \emph{BstBI},
#'       \emph{EcoT22I}, \emph{FspI}, \emph{MluI-HF}, \emph{NruI-HF}, \emph{PstI},
#'       and \emph{SbfI}.}
#'   \item{sites}{Enzyme binding site sequences. This vector is of the binding sites,
#'       5' then 3', for each unique site that the enzyme can bind to.}
#' }
#' @source \url{http://www.neb.com}
"binding_sites"


#' Table of evolutionary rates.
#'
#' From...
#'
#' Sung, W., M. S. Ackerman, M. M. Dillon, T. G. Platt, C. Fuqua, V. S. Cooper, and
#' M. Lynch. 2016. Evolution of the insertion-deletion mutation rate across the
#' tree of life. \emph{G3: Genes|Genomes|Genetics} \strong{6}:2583–2591.
#'
#'
#' @format A data frame with 15 rows and 4 variables:
#' \describe{
#'   \item{domain}{Either \code{Bacteria} or \code{Eukarya} for what type of organism
#'       the species is.}
#'   \item{species}{Species name.}
#'   \item{indels}{Rate of insertions and deletions (events per site per generation).}
#'   \item{subs}{Base-substitution mutation rate (events per site per generation).}
#' }
#'
#' @source \url{http://dx.doi.org/10.1534/g3.116.030890}
#'
"evo_rates"




