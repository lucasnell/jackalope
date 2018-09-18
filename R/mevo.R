
#' Make substitution rate matrix (Q) and vector of equilibrium frequencies (pis).
#'
#' @inheritParams make_mevo
#'
#' @noRd
#'
#'
substitutions <- function(sub) {

    err_msg <- paste("\nNot all required parameters provided for the given",
                     "substitution model (\"%s\").",
                     "See `vignette(\"sub-model\")` for each model's required",
                     "parameter(s).")

    if (sub$model == "TN93") {
        if (any(! c("pi_tcag", "alpha_1", "alpha_2", "beta") %in% names(sub))) {
            stop(sprintf(err_msg, sub$model), call. = FALSE)
        }
        Q <- TN93_rate_matrix(sub$pi_tcag, sub$alpha_1, sub$alpha_2,
                              sub$beta)
        pi_tcag <- sub$pi_tcag
    } else if (sub$model == "JC69") {
        if (any(! c("lambda") %in% names(sub))) {
            stop(sprintf(err_msg, sub$model), call. = FALSE)
        }
        Q <- JC69_rate_matrix(sub$lambda)
        pi_tcag <- rep(0.25, 4)
    } else if (sub$model == "K80") {
        if (any(! c("alpha", "beta") %in% names(sub))) {
            stop(sprintf(err_msg, sub$model), call. = FALSE)
        }
        Q <- K80_rate_matrix(sub$alpha, sub$beta);
        pi_tcag <- rep(0.25, 4)
    } else if (sub$model == "F81") {
        if (any(! c("pi_tcag") %in% names(sub))) {
            stop(sprintf(err_msg, sub$model), call. = FALSE)
        }
        Q <- F81_rate_matrix(sub$pi_tcag);
        pi_tcag <- sub$pi_tcag
    } else if (sub$model == "HKY85") {
        if (any(! c("pi_tcag", "alpha", "beta") %in% names(sub))) {
            stop(sprintf(err_msg, sub$model), call. = FALSE)
        }
        Q <- HKY85_rate_matrix(sub$pi_tcag, sub$alpha, sub$beta)
        pi_tcag <- sub$pi_tcag
    } else if (sub$model == "F84") {
        if (any(! c("pi_tcag", "beta", "kappa") %in% names(sub))) {
            stop(sprintf(err_msg, sub$model), call. = FALSE)
        }
        Q <- F84_rate_matrix(sub$pi_tcag, sub$beta, sub$kappa)
        pi_tcag <- sub$pi_tcag
    } else if (sub$model == "GTR") {
        if (any(! c("pi_tcag", "abcdef") %in% names(sub))) {
            stop(sprintf(err_msg, sub$model), call. = FALSE)
        }
        Q <- GTR_rate_matrix(sub$pi_tcag, sub$abcdef)
        pi_tcag <- sub$pi_tcag
    } else if (sub$model == "UNREST") {
        if (any(! c("Q") %in% names(sub))) {
            stop(sprintf(err_msg, sub$model), call. = FALSE)
        }
        q_pi_list <- UNREST_rate_matrix(sub$Q)
        Q <- q_pi_list$Q
        pi_tcag <- q_pi_list$pi_tcag
    } else stop("\nInvalid `sub$model` argument to `make_mevo`.", call. = FALSE)

    return(list(Q = Q, pi_tcag = pi_tcag))
}




#' Process information for indels
#'
#' @param indel List of indel parameters.
#' @param which_type String of `"sub"` or `"del"` to print for error messages.
#'
#' @noRd
#'
indels <- function(indel) {

    which_type <- deparse(substitute(indel))

    if (is.null(indel)) {
        rates <- numeric(0)
    } else {
        err_msg <- paste("\nWhen specifying", which_type, "in `make_mevo`, ")
        if (is.null(indel$rate)) {
            stop(err_msg, "you must always provide a rate.",
                 call. = TRUE)
        } else if (!single_number(indel$rate, 0)) {
            stop(err_msg, "the rate must be a single number >= 0.",
                 call. = TRUE)
        }
        if (indel$rate == 0) return(numeric(0))
        names_ <- sort(names(indel))
        if (identical(names_, sort(c("rate", "max_length")))) {
            if (!single_whole_number(indel$max_length, 1)) {
                stop(err_msg, "the max length must be a single whole number >= 1.",
                     call. = TRUE)
            }
            rel_rates <- exp(-1 * 1:(indel$max_length))
        } else if (identical(names_, sort(c("rate", "max_length", "a")))) {
            if (!single_whole_number(indel$max_length, 1)) {
                stop(err_msg, "the max length must be a single whole number >= 1.",
                     call. = TRUE)
            }
            L <- 1:(indel$max_length)
            rel_rates <- {(L * indel$max_length) / (indel$max_length - L + 1)}^(-indel$a)
        } else if (identical(names_, sort(c("rate", "rel_rates")))) {
            rel_rates <- indel$rel_rates
        } else {
            stop(err_msg, "it must contain names that coincide with one of the methods ",
                 "in the \"Indels\" section within `?make_mevo`.",
                 "Note that extra names return an error.",
                 call. = FALSE)
        }

        # So relative rates sum to 1:
        rel_rates <- rel_rates / sum(rel_rates)

        # Absolute rates:
        rates <- rel_rates * indel$rate

    }

    return(rates)
}




#' Process information for variation in mutation rates among sites.
#'
#' @inheritParams make_mevo
#' @param seq_sizes Vector of reference-genome sequence sizes.
#'
#'
#' @noRd
#'
site_variability <- function(site_var, seq_sizes) {

    if (!is.null(site_var)) {

        if (all(c("shape", "region_size", "mats") %in% names(site_var))) {
            stop("\nThe `site_var` argument to `make_mevo` ",
                 "must be a named list with names of either ",
                 "(a) both \"shape\" and \"region_size\" or ",
                 "(b) just \"mats\". ",
                 "Providing all three names is not permitted.",
                 call. = FALSE)
        }

        if (all(c("shape", "region_size") %in% names(site_var))) {
            gamma_mats <- make_gamma_mats(seq_sizes,
                                          gamma_size_ = site_var$region_size,
                                          shape = site_var$shape)
        } else if ("mats" %in% names(site_var)) {

            err_msg <- paste("\nThe `mats` field inside the `site_var`",
                             "argument to the `make_mevo` function needs to",
                             "be a list of matrices.")
            if (!inherits(site_var$mats, "list")) {
                stop(err_msg, call. = FALSE)
            } else if (!all(sapply(site_var$mats, inherits,
                                   what = "matrix"))) {
                stop(err_msg, call. = FALSE)
            }

            # Check matrices for proper end points and # columns:
            check_gamma_mats(site_var$mats, seq_sizes)

            gamma_mats <- site_var$mats

        } else {
            stop("\nThe `site_var` argument to `make_mevo` ",
                 "must be a named list with names of either ",
                 "(a) both \"shape\" and \"region_size\" or ",
                 "(b) just \"mats\".",
                 call. = FALSE)
        }
    } else {
        # This results in no variability among sites:
        gamma_mats <- make_gamma_mats(seq_sizes, gamma_size_ = 0, shape = 1)
    }
    return(gamma_mats)
}




#' Make a `mevo` object to store information needed for molecular evolution simulation.
#'
#'
#'
#' @section Indels:
#' Both insertions and deletions require the `rate` parameter, which specifies
#' the overall insertion/deletion rate among all lengths.
#' The `rate` parameter is ultimately combined with a vector of relative rates among
#' the different lengths of insertions/deletions from 1 to the maximum
#' possible length.
#' There are three different ways to specify/generate relative-rate values.
#' \enumerate{
#'     \item Assume that rates are proportional to `exp(-L)` for indel length
#'         `L` from 1 to the maximum length (Albers et al. 2011).
#'         This method will be used if the following arguments are provided:
#'         \itemize{
#'             \item `rate`
#'             \item `max_length`
#'         }
#'     \item Generate relative rates from a Lavalette distribution
#'         (Fletcher and Yang 2009), where the rate for length `L` is proportional to
#'         `{L * max_size / (max_size - L + 1)}^(-a)`.
#'         This method will be used if the following arguments are provided:
#'         \itemize{
#'             \item `rate`
#'             \item `max_length`
#'             \item `a`
#'         }
#'     \item Directly specify values by providing a numeric vector of relative
#'         rates for each insertion/deletion length from 1 to the maximum length.
#'         This method will be used if the following arguments are provided:
#'         \itemize{
#'             \item `rate`
#'             \item `rel_rates`
#'         }
#' }
#' Note that providing extra parameter(s) will return an error.
#'
#'
#' @section Site variation:
#' A site's deviance from the average mutation rate is determined by its
#' "gamma distance".
#' A site's overall mutation rate is the mutation rate for that nucleotide
#' (substitution + indel) multiplied by the site's gamma distance.
#' There are two options for specifying gamma distances:
#' \enumerate{
#'     \item Generate gamma distances from a Gamma distribution.
#'         This option requires the following arguments:
#'         \describe{
#'             \item{`shape`}{Shape parameter for the Gamma distribution,
#'                 where the variance of the distribution is `1 / shape`.
#'                 The mean is fixed to 1.}
#'             \item{`region_size`}{Size of regions to break the genome into,
#'                 where all sites within a region have the same gamma distance.}
#'         }
#'     \item Manually input matrices that specify the gamma distance and end points
#'         for regions each gamma distances refers to.
#'         This option requires the following argument:
#'         \describe{
#'             \item{`mats`}{List of matrices, one for each sequence in the genome.
#'                 Each matrix should have two columns.
#'                 The first should contain the end points for each region.
#'                 The second should contain the gamma distances for each region.
#'                 Note that if gamma distances don't have a mean (weighted by
#'                 sequence length for each gamma-distance value) equal to 1,
#'                 you're essentially changing the overall mutation rate.}
#'         }
#' }
#'
#'
#' @param reference A \code{ref_genome} object from which you will generate variants.
#' @param sub A list containing the parameters for substitutions.
#'     The `model` field is always required within this list, and this specifies
#'     the model used for substitutions.
#'     It takes one of the following options:
#'     `"JC69"`, `"K80"`, `"F81"`, `"HKY85"`, `"TN93"`, `"F84"`, `"GTR"`, or `"UNREST"`.
#'     See `vignette("sub-models")` for more information on these models and
#'     their required parameters.
#' @param ins A list containing the parameters for insertions.
#'     Passing `NULL` to this argument results in no insertions.
#'     See "Indels" section for more information.
#'     Defaults to `NULL`.
#' @param del A list containing the parameters for deletions.
#'     Passing `NULL` to this argument results in no deletions.
#'     See "Indels" section for more information.
#'     Defaults to `NULL`.
#' @param site_var List of parameters for generating variability in mutation
#'     rates among sites (for both substitutions and indels).
#'     Passing `NULL` to this argument results in no variability among sites.
#'     See "Site variation" section for more information.
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
#' @references
#' Albers, C. A., G. Lunter, D. G. MacArthur, G. McVean, W. H. Ouwehand, and R. Durbin.
#' 2011. Dindel: accurate indel calls from short-read data. Genome Research 21:961–973.
#'
#' Fletcher, W., and Z. Yang. 2009. INDELible: a flexible simulator of
#' biological sequence evolution. Molecular Biology and Evolution 26:1879–1888.
#'
#' @examples
#' ref <- create_genome(10, 1000)
#' mevo <- make_mevo(ref, list(model = "JC69", lambda = 1))
#'
make_mevo <- function(reference,
                      sub,
                      ins = NULL,
                      del = NULL,
                      site_var = NULL,
                      chunk_size = 100) {

    if (!inherits(reference, "ref_genome")) {
        stop("\nIn `make_mevo`, the `reference` argument must be a ref_genome object.",
             call. = FALSE)
    }

    # -------+
    # Process substitution info:
    # -------+
    if (is.null(sub$model)) {
        stop("\nYou must always specify a `model` field to the `sub` argument in ",
             "`make_mevo` (i.e., `is.null(sub$model)` should never be `TRUE`).",
             call. = FALSE)
    }
    sub$model <- match.arg(sub$model, c("JC69", "K80", "F81", "HKY85", "TN93", "F84",
                                        "GTR", "UNREST"))
    sub_info <- substitutions(sub)


    # -------+
    # Process indel info:
    # -------+
    insertion_rates <- indels(ins)
    deletion_rates <- indels(del)


    # -------+
    # Process info for mutation-rate variability among sites:
    # -------+
    gamma_mats <- site_variability(site_var, seq_sizes = reference$sizes())


    # -------+
    # Make final output object
    # -------+
    out <- mevo$new(sub_info,
                    insertion_rates,
                    deletion_rates,
                    gamma_mats,
                    chunk_size)

    return(out)
}

