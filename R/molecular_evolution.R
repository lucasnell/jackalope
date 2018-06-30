
#' Make substitution rate matrix (Q) and vector of equilibrium frequencies (pis).
#'
#' @inheritParams make_sampler
#' @param xi Overall indel rate.
#'
#' @noRd
#'
#'
make_Q_pis <- function(sub_params, model, xi) {

    err_msg <- paste("\nNot all required names provided in `sub_params`.",
                     "See `?make_samplers` for what to provide.")

    if (model == "TN93") {
        if (any(! c("pi_tcag", "alpha_1", "alpha_2", "beta") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        Q <- TN93_rate_matrix(sub_params$pi_tcag, sub_params$alpha_1, sub_params$alpha_2,
                              sub_params$beta, xi)
        pi_tcag <- sub_params$pi_tcag
    } else if (model == "JC69") {
        if (any(! c("lambda") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        Q <- JC69_rate_matrix(sub_params$lambda, xi)
        pi_tcag <- sub_params$pi_tcag
    } else if (model == "K80") {
        if (any(! c("alpha", "beta") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        Q <- K80_rate_matrix(sub_params$alpha, sub_params$beta, xi);
        pi_tcag <- sub_params$pi_tcag
    } else if (model == "F81") {
        if (any(! c("pi_tcag") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        Q <- F81_rate_matrix(sub_params$pi_tcag, xi);
        pi_tcag <- sub_params$pi_tcag
    } else if (model == "HKY85") {
        if (any(! c("pi_tcag", "alpha", "beta") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        Q <- HKY85_rate_matrix(sub_params$pi_tcag, sub_params$alpha, sub_params$beta,
                               xi)
        pi_tcag <- sub_params$pi_tcag
    } else if (model == "F84") {
        if (any(! c("pi_tcag", "beta", "kappa") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        Q <- F84_rate_matrix(sub_params$pi_tcag, sub_params$beta, sub_params$kappa,
                             xi)
        pi_tcag <- sub_params$pi_tcag
    } else if (model == "GTR") {
        if (any(! c("pi_tcag", "abcdef") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        Q <- GTR_rate_matrix(sub_params$pi_tcag, sub_params$abcdef, xi)
        pi_tcag <- sub_params$pi_tcag
    } else if (model == "UNREST") {
        if (any(! c("Q") %in% names(sub_params))) {
            stop(err_msg, call. = FALSE)
        }
        q_pi_list <- UNREST_rate_matrix(sub_params$Q, xi)
        Q <- q_pi_list$Q
        pi_tcag <- q_pi_list$pi_tcag
    } else stop("\nInput model to `make_sampler` not available.", call. = FALSE)

    return(list(Q = Q, pi_tcag = pi_tcag))
}

#' Create mutation sampler objects.
#'
#' @param sub_params A list containing the necessary parameters for the specified
#'     substitution model.
#' @param indel_params A list containing the necessary parameters for indels.
#' It requires the following parameters:
#' \describe{
#' \item{`xi`}{Overall indel rate.}
#' \item{`psi`}{Proportion of insertions to deletions.}
#' \item{`rel_insertion_rates`}{Relative insertion rates.}
#' \item{`rel_deletion_rates`}{Relative deletion rates.}
#' }
#' @param chunk_size The size of chunks to first sample uniformly before doing
#'     weighted sampling by rates.
#'     Uniformly sampling before doing weighted sampling dramatically speeds up
#'     the mutation process (especially for very long sequences) and has little
#'     effect on the sampling probabilities.
#'     Set this to `0` to not uniformly sample first.
#'     Defaults to `100`.
#' @param model Character indicating which substitution mutation model to use.
#' Takes one of the following options, with the required parameters in parentheses:
#'
#' - `"TN93"` (`pi_tcag`, `alpha_1`, `alpha_2`, `beta`)
#' - `"JC69"` (`lambda`)
#' - `"K80"` (`alpha`, `beta`)
#' - `"F81"` (`pi_tcag`)
#' - `"HKY85"` (`pi_tcag`, `alpha`, `beta`)
#' - `"F84"` (`pi_tcag`, `beta`, `kappa`)
#' - `"GTR"` (`pi_tcag`, `abcdef`)
#' - `"UNREST"` (`Q`)
#'
#' Defaults to `"TN93"`.
#'
#'
#' @noRd
#'
make_sampler_ptr <- function(sub_params,
                             indel_params,
                             model = c("TN93", "JC69", "K80", "F81", "HKY85",
                                       "F84", "GTR", "UNREST"),
                             chunk_size = 100) {

    if (any(! c("xi", "psi", "rel_insertion_rates", "rel_deletion_rates") %in%
            names(indel_params))) {
        stop("\nNot all required names provided in `indel_params`. ",
             "See `?make_samplers` for what to provide.", call. = FALSE)
    }

    model <- match.arg(model)

    Q_pis <- make_Q_pis(sub_params, model, indel_params$xi)

    if (chunk_size <= 0) {
        sampler_ptr <- make_mutation_sampler_base(Q_pis$Q,
                                                  indel_params$xi,
                                                  indel_params$psi,
                                                  Q_pis$pi_tcag,
                                                  indel_params$rel_insertion_rates,
                                                  indel_params$rel_deletion_rates)
    } else {
        sampler_ptr <- make_mutation_sampler_chunk_base(Q_pis$Q,
                                                        indel_params$xi,
                                                        indel_params$psi,
                                                        Q_pis$pi_tcag,
                                                        indel_params$rel_insertion_rates,
                                                        indel_params$rel_deletion_rates,
                                                        chunk_size)
    }

    return(sampler_ptr)

}

