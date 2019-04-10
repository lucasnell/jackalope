
#' Table of evolutionary rates.
#'
#' From Table 1 in Sung et al. (2016).
#'
#'
#' @format A data frame with 15 rows and 8 variables:
#' \describe{
#'   \item{`domain`}{Either \code{Bacteria} or \code{Eukarya} for what type of organism
#'       the species is.}
#'   \item{`species`}{Species name.}
#'   \item{`Ge`}{Effective genome size using only coding DNA.}
#'   \item{`Gc_Gnc`}{Effective genome size using coding DNA and non-coding DNA that is
#'       under purifying selection.}
#'   \item{`indels`}{Rate of insertions and deletions (\eqn{10^{-10}}{10^-10}
#'       events per site per generation).}
#'   \item{`subs`}{Base-substitution mutation rate (\eqn{10^{-10}}{10^-10}
#'       events per site per generation).}
#'   \item{`Ne`}{Effective population size (\eqn{\times 10^{6}}{* 10^6}).}
#'   \item{`theta_s`}{Population mutation rate estimated using \eqn{\theta_s}{theta_s}.}
#'   \item{`pi_s`}{Population mutation rate estimated using \eqn{\pi_s}{pi_s}.}
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


