
#' Make substitution rate matrix (Q) and vector of equilibrium frequencies (pis).
#'
#' @inheritParams make_mevo
#' @param xi Overall indel rate.
#'
#' @noRd
#'
#'
make_Q_pis <- function(sub_params, sub_model, xi) {

    err_msg <- paste("\nNot all required names provided in `sub_params`.",
                     "See `?create_variants` for what to provide for each",
                     "possible value of the `sub_model` argument.")

    if (sub_model == "TN93") {
        if (any(! c("pi_tcag", "alpha_1", "alpha_2", "beta") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        Q <- TN93_rate_matrix(sub_params$pi_tcag, sub_params$alpha_1, sub_params$alpha_2,
                              sub_params$beta, xi)
        pi_tcag <- sub_params$pi_tcag
    } else if (sub_model == "JC69") {
        if (any(! c("lambda") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        Q <- JC69_rate_matrix(sub_params$lambda, xi)
        pi_tcag <- sub_params$pi_tcag
    } else if (sub_model == "K80") {
        if (any(! c("alpha", "beta") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        Q <- K80_rate_matrix(sub_params$alpha, sub_params$beta, xi);
        pi_tcag <- sub_params$pi_tcag
    } else if (sub_model == "F81") {
        if (any(! c("pi_tcag") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        Q <- F81_rate_matrix(sub_params$pi_tcag, xi);
        pi_tcag <- sub_params$pi_tcag
    } else if (sub_model == "HKY85") {
        if (any(! c("pi_tcag", "alpha", "beta") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        Q <- HKY85_rate_matrix(sub_params$pi_tcag, sub_params$alpha, sub_params$beta,
                               xi)
        pi_tcag <- sub_params$pi_tcag
    } else if (sub_model == "F84") {
        if (any(! c("pi_tcag", "beta", "kappa") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        Q <- F84_rate_matrix(sub_params$pi_tcag, sub_params$beta, sub_params$kappa,
                             xi)
        pi_tcag <- sub_params$pi_tcag
    } else if (sub_model == "GTR") {
        if (any(! c("pi_tcag", "abcdef") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        Q <- GTR_rate_matrix(sub_params$pi_tcag, sub_params$abcdef, xi)
        pi_tcag <- sub_params$pi_tcag
    } else if (sub_model == "UNREST") {
        if (any(! c("Q") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        q_pi_list <- UNREST_rate_matrix(sub_params$Q, xi)
        Q <- q_pi_list$Q
        pi_tcag <- q_pi_list$pi_tcag
    } else stop("\nInvalid `sub_model` argument to `create_variants`.", call. = FALSE)

    return(list(Q = Q, pi_tcag = pi_tcag))
}











#' Make a `mevo` object to store information needed for molecular evolution simulation.
#'
#'
#' @param sub_model Character indicating which substitution mutation model to use.
#'     It takes one of the following options, with the optional parameters in parentheses
#'     for each:
#'     \describe{
#'         \item{`"TN93"`}{ `pi_tcag`, `alpha_1`, `alpha_2`, `beta` }
#'         \item{`"JC69"`}{ `lambda` }
#'         \item{`"K80"`}{ `alpha`, `beta` }
#'         \item{`"F81"`}{ `pi_tcag` }
#'         \item{`"HKY85"`}{ `pi_tcag`, `alpha`, `beta` }
#'         \item{`"F84"`}{ `pi_tcag`, `beta`, `kappa` }
#'         \item{`"GTR"`}{ `pi_tcag`, `abcdef` }
#'         \item{`"UNREST"`}{ `Q` }
#'     }
#'     See `?vignette("sub-models")` for more information on these model parameters
#'     and their default values.
#'     Defaults to `"TN93"`.
#' @param sub_params A list containing the parameters for the specified
#'     substitution model.
#'     See `?vignette("sub-models")` for more information on these model parameters
#'     and their default values.
#'     Defaults to `NULL`, which causes it to use all default parameters.
#' @param indel_params A list containing the parameters for indels.
#'     The following parameters are allowed:
#'     \describe{
#'         \item{`xi`}{Overall indel rate.}
#'         \item{`psi`}{Proportion of insertions to deletions.}
#'         \item{`rel_insertion_rates`}{Relative insertion rates.}
#'         \item{`rel_deletion_rates`}{Relative deletion rates.}
#'     }
#'     Defaults to `NULL`, which causes it to use all default parameters. See
#'     Details for default parameter values.
#' @param site_var_params List of parameters for generating variability in mutation
#'     rates among sites (for both substitutions and indels).
#'     A site's deviance from the average mutation rate is determined by its
#'     "gamma distance".
#'     A site's overall mutation rate is the mutation rate for that nucleotide
#'     (substitution + indel) multiplied by the site's gamma distance.
#'     There are two options for specifying gamma distances:
#'     \enumerate{
#'         \item Generate gamma distances from a Gamma distribution.
#'             This option requires the following arguments:
#'             \describe{
#'                 \item{`shape`}{Shape parameter for the Gamma distribution,
#'                     where the variance of the distribution is `1 / shape`.
#'                     The mean is fixed to 1.}
#'                 \item{`region_size`}{Size of regions where each site within that
#'                     region has the same gamma distance.}
#'             }
#'         \item Manually input matrices that specify the gamma distance and end points
#'             for regions each gamma distances refers to.
#'             This option requires the following argument:
#'             \describe{
#'                 \item{`mats`}{List of matrices, one for each sequence in the genome.
#'                     Each matrix should have two columns.
#'                     The first should contain the end points for each region.
#'                     The second should contain the gamma distances for each region.
#'                     Note that if gamma distances don't have a mean (weighted by
#'                     sequence length for each gamma-distance value) equal to 1,
#'                     you're essentially changing the overall mutation rate.}
#'             }
#'     }
#'     Passing `NULL` to this argument results in no variability among sites.
#'     Defaults to `NULL`.
#' @param chunk_size The size of "chunks" of sequences to first sample uniformly
#'     before doing weighted sampling by rates for each sequence location.
#'     Uniformly sampling before doing weighted sampling dramatically speeds up
#'     the mutation process (especially for very long sequences) and has little
#'     effect on the sampling probabilities.
#'     Higher values will more closely resemble sampling without the uniform-sampling
#'     step, but will be slower.
#'     Set this to `0` to not uniformly sample first.
#'     From testing on a chromosome of length `1e6`, a `chunk_size` value of `100`
#'     offers a ~10x speed increase and doesn't differ significantly from sampling
#'     without the uniform-sampling step.
#'     Defaults to `100`.
#'
#' @return An object of class \code{\link{mevo}}.
#'
#' @export
#'
#' @examples
#'
make_mevo <- function(reference,
                      sub_model = "TN93",
                      sub_params = NULL,
                      indel_params = NULL,
                      site_var_params = NULL,
                      chunk_size = 100) {

    if (!inherits(reference, "ref_genome")) {
        stop("\nIn `make_mevo`, the `reference` argument must be a ref_genome object.",
             call. = FALSE)
    }

    sub_model <- match.arg(sub_model, c("TN93", "JC69", "K80", "F81", "HKY85",
                                        "F84", "GTR", "UNREST"))

    # ================================*
    # ================================*

    # ** LEFT OFF ----
    # Manage creation of matrices and indel vectors in a more elegant way.

    # ================================*
    # ================================*


    # * start temp check ----
    if (is.null(sub_params) | is.null(indel_params)) {
        stop("\nIn make_mevo, defaults are not yet programmed for sub_params and ",
             "indel_params, so you cannot pass NULL for either of these arguments.",
             call. = FALSE)
    }
    # _ end temp check ----



    if (any(! c("xi", "psi", "rel_insertion_rates", "rel_deletion_rates") %in%
            names(indel_params))) {
        stop("\nNot all required names provided in `indel_params`. ",
             "See `?make_mevo` for what to provide.", call. = FALSE)
    }

    Q_pis <- make_Q_pis(sub_params, sub_model, indel_params$xi)

    # * start temp assign ----
    # Below assigns these values to objects in global namespace.
    # This is not a permanent solution.
    xi <- indel_params$xi
    psi <- indel_params$psi
    rel_insertion_rates <- indel_params$rel_insertion_rates
    rel_deletion_rates <- indel_params$rel_deletion_rates
    # _ end temp assign ----



    # -------+
    # Make Gamma matrices (for mutation-rate variability among sites):
    # -------+
    seq_sizes <- reference$sizes
    if (!is.null(site_var_params)) {
        if (all(c("shape", "region_size") %in% names(site_var_params))) {
            gamma_mats <- make_gamma_mats(seq_sizes,
                                          gamma_size_ = site_var_params$region_size,
                                          shape = site_var_params$shape)
        } else if ("mats" %in% names(site_var_params)) {

            err_msg <- paste("\nThe `mats` field inside the `site_var_params`",
                             "argument to the `create_variants` function needs to",
                             "be a list of matrices.")
            if (!inherits(site_var_params$mats, "list")) {
                stop(err_msg, call. = FALSE)
            } else if (!all(sapply(site_var_params$mats, inherits,
                                   what = "matrix"))) {
                stop(err_msg, call. = FALSE)
            }

            # Check matrices for proper end points and # columns:
            check_gamma_mats(site_var_params$mats, seq_sizes)

            gamma_mats <- site_var_params$mats

        } else {
            stop("\nThe `site_var_params` argument to `create_variants` ",
                 "must be a named list containing \"shape\" and \"region_size\" ",
                 "or \"mats\" as names.",
                 call. = FALSE)
        }
    } else {
        # This results in no variability among sites:
        gamma_mats <- make_gamma_mats(seq_sizes, gamma_size_ = 0, shape = 1)
    }

    out <- mevo$new(Q_pis$Q,
                    xi,
                    psi,
                    Q_pis$pi_tcag,
                    rel_insertion_rates,
                    rel_deletion_rates,
                    gamma_mats,
                    chunk_size)

    return(out)
}

