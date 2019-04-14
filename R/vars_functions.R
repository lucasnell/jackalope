

#' Organize higher-level information for creating variants.
#'
#'
#' The following functions organize information that gets passed to `create_variants`
#' to generate variants from a reference genome.
#' Each function represents a method of generation and starts with `"vars_"`.
#' The first three are phylogenomic methods, and all functions but `vars_vcf`
#' will use molecular evolution information when passed to `create_variants`.
#'
#' \describe{
#'     \item{\code{\link{vars_theta}}}{Uses an estimate for theta, the population-scaled
#'         mutation rate, and a desired number of variants.}
#'     \item{\code{\link{vars_phylo}}}{Uses phylogenetic tree(s) from `phylo`
#'         object(s) or NEWICK file(s), one tree per sequence or one for all sequences.}
#'     \item{\code{\link{vars_gtrees}}}{Uses gene trees, either in the form of
#'         an object from the `scrm` or `coala` package or
#'         a file containing output in the style of the `ms` program.}
#'     \item{\code{\link{vars_ssites}}}{Uses matrices of segregating sites,
#'         either in the form of
#'         `scrm` or `coala` coalescent-simulator object(s), or
#'         (2) a `ms`-style output file.}
#'     \item{\code{\link{vars_vcf}}}{Uses a variant call format (VCF) file that
#'         directly specifies variants.
#'         This method does not work if the `vcfR` package isn't installed.}
#' }
#'
#'
#' @seealso \code{\link{create_variants}}
#'
#' @name vars_functions
#'
NULL
