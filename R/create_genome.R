

#' Create a reference genome.
#'
#' Random chromosomes are generated to create a new \code{ref_genome} object.
#' Note that this function will never generate empty chromosomes.
#'
#' @param n_chroms Number of chromosomes.
#' @param len_mean Mean for the gamma distribution of chromosome sizes.
#' @param len_sd Standard deviation for the gamma distribution of chromosome sizes.
#'     If set to \code{<= 0}, all chromosomes will be the same length. Defaults to \code{0}.
#' @param pi_tcag Vector of length 4 containing the nucleotide equilibrium frequencies
#'     for "T", "C", "A", and "G", respectively. Defaults to `rep(0.25, 4)`.
#' @param n_threads Number of threads to use for parallel processing. This argument is
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
create_genome <- function(n_chroms,
                          len_mean,
                          len_sd = 0,
                          pi_tcag = rep(0.25, 4),
                          n_threads = 1) {


    if (!single_integer(n_chroms, .min = 1)) {
        err_msg("create_genome", "n_chroms", "a single integer >= 1")
    }
    if (!single_number(len_mean, .min = 1)) {
        err_msg("create_genome", "len_mean", "a single number >= 1")
    }
    if (!single_number(len_sd, .min = 0)) {
        err_msg("create_genome", "len_sd", "a single number >= 0")
    }

    if (!is_type(pi_tcag, "numeric", 4) || any(pi_tcag < 0) || all(pi_tcag == 0)) {
        err_msg("create_genome", "pi_tcag", "a numeric vector of length 4,",
                "where no number can be < 0 and at least one must be > 0")
    }
    if (!single_integer(n_threads, .min = 1)) {
        err_msg("create_genome", "n_threads", "a single integer >= 1")
    }

    ptr <- create_genome_cpp(n_chroms, len_mean, len_sd, pi_tcag, n_threads)

    ref_obj <- ref_genome$new(ptr)

    return(ref_obj)
}
