

#' Create a reference genome.
#'
#' Random sequences are generated to create a new \code{ref_genome} object.
#' Note that this function will never generate empty sequences.
#'
#' @param n_seqs Number of sequences.
#' @param len_mean Mean for the gamma distribution of sequence sizes.
#' @param len_sd Standard deviation for the gamma distribution of sequence sizes.
#'     If set to \code{<= 0}, all sequences will be the same length. Defaults to \code{0}.
#' @param pi_tcag Vector of length 4 containing the nucleotide equilibrium frequencies
#'     for "T", "C", "A", and "G", respectively. Defaults to `rep(0.25, 4)`.
#' @param n_cores Number of cores to use for parallel processing. This argument is
#'     ignored if OpenMP is not enabled. Defaults to \code{1}.
#'
#'
#' @return A \code{\link{ref_genome}} object.
#'
#' @export
#'
#' @examples
#'
#' genome <- create_genome(10, 100e3, 100, pi_tcag = c(0.1, 0.2, 0.3, 0.4))
#'
create_genome <- function(n_seqs,
                          len_mean,
                          len_sd = 0,
                          pi_tcag = rep(0.25, 4),
                          n_cores = 1) {


    if (!single_whole_number(n_seqs, .min = 1)) {
        stop("\nThe n_seqs argument supplied to create_genome is not a single",
             "whole number greater than 0.",
             call. = FALSE)
    }
    if (!single_number(len_mean, .min = 1)) {
        stop("\nThe len_mean argument supplied to create_genome is not a single",
             "number greater than or equal to 1.",
             call. = FALSE)
    }
    if (!single_number(len_sd, .min = 0)) {
        stop("\nThe len_sd argument supplied to create_genome is not a single",
             "number >= 0.",
             call. = FALSE)
    }

    if (!is.numeric(pi_tcag) | length(pi_tcag) != 4 | any(pi_tcag < 0)) {
        stop("\nThe pi_tcag argument supplied to create_genome is not a numeric",
             "vector of exactly 4 non-negative numbers.",
             call. = FALSE)
    }
    if (!single_whole_number(n_cores, .min = 1)) {
        stop("\nThe n_cores argument supplied to create_genome is not a single",
             "whole number greater than or equal to 1.",
             call. = FALSE)
    }

    ptr <- create_genome_(n_seqs, len_mean, len_sd, pi_tcag, n_cores)

    ref_obj <- ref_genome$new(ptr)

    return(ref_obj)
}
