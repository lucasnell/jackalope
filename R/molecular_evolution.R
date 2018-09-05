
#' Make substitution rate matrix (Q) and vector of equilibrium frequencies (pis).
#'
#' @inheritParams make_sampler_ptr
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




#' Create mutation sampler objects.
#'
#' @inheritParams create_variants
#'
#' @noRd
#'
make_sampler_ptr <- function(sub_params,
                             indel_params,
                             sub_model = c("TN93", "JC69", "K80", "F81", "HKY85",
                                       "F84", "GTR", "UNREST"),
                             chunk_size = 100) {

    if (any(! c("xi", "psi", "rel_insertion_rates", "rel_deletion_rates") %in%
            names(indel_params))) {
        stop("\nNot all required names provided in `indel_params`. ",
             "See `?create_variants` for what to provide.", call. = FALSE)
    }

    Q_pis <- make_Q_pis(sub_params, sub_model, indel_params$xi)

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

