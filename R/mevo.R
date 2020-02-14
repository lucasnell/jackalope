
#' ========================================================================`
#'
#' This file stores functions for molecular-evolution info.
#' This includes those for substitutions and indels.
#'
#' ========================================================================`


# SUBSTITUTIONS -----


# ... sub_info class -------



#' An R6 class representing information for a substitution model.
#'
#'
#' This class should NEVER be created using `sub_info$new`.
#' Only use one of the functions in `?sub_models`.
#' That's why I'm not exporting it.
#'
#' @noRd
#'
#'
#' @importFrom R6 R6Class
#'
sub_info <- R6Class(

    "sub_info",

    public = list(

        initialize = function(info_list) {

            extra_msg <- paste(" Please only create these objects using one of the",
                               "functions in ?sub_models, NOT using sub_info$new().")
            if (!inherits(info_list, "list")) {
                stop("\nWhen initializing a sub_info object, you need to use ",
                     "a list.", extra_msg, call. = FALSE)
            }

            for (x in c("Q", "pi_tcag", "U", "Ui", "L", "gammas",
                        "invariant", "model")) {
                if (is.null(info_list[[x]])) {
                    stop("\nWhen initializing a sub_info object, the input list ",
                         "must contain the field ", x, ".", extra_msg, call. = FALSE)
                }
            }

            private$r_Q <- info_list$Q
            private$r_pi_tcag <- info_list$pi_tcag
            private$r_U <- info_list$U
            private$r_Ui <- info_list$Ui
            private$r_L <- info_list$L
            private$r_gammas <- info_list$gammas
            private$r_invariant <- info_list$invariant
            private$r_model <- info_list$model

            # Revert these back to list (from arma::field which adds dims):
            dim(private$r_Q) <- NULL
            dim(private$r_U) <- NULL
            dim(private$r_Ui) <- NULL
            dim(private$r_L) <- NULL
        },

        print = function(...) {

            digits <- max(3, getOption("digits") - 3)

            cat(sprintf("< Substitution model %s >\n", private$r_model))
            fmt <- paste0("%.", digits, "f")

            cat("# Equilibrium densities:\n")
            cat("  ", sprintf(fmt, private$r_pi_tcag), "\n")

            cat("# Substitution rate matrix:\n")
            prmatrix(private$r_Q[[1]] / private$r_gammas[1], digits = digits,
                     rowlab = paste("  ", c("T", "C", "A", "G")),
                     collab = c("T", "C", "A", "G"))

            if (length(private$r_gammas) == 1 & private$r_invariant == 0) {
                cat("# No among-site variability\n")
            } else {
                if (length(private$r_gammas) > 1) {
                    cat(sprintf("# Discrete Gamma classes: %i\n", length(private$r_gammas)))
                } else {
                    cat("# No continuous variability among sites\n")
                }
                if (private$r_invariant > 0) {
                    cat(sprintf(paste("# Proportion of invariant sites:", fmt, "\n"),
                                private$r_invariant))
                } else {
                    cat("# No invariant sites\n")
                }
            }

            invisible(self)

        },

        Q = function() return(private$r_Q),
        pi_tcag = function() return(private$r_pi_tcag),
        U = function() return(private$r_U),
        Ui = function() return(private$r_Ui),
        L = function() return(private$r_L),
        gammas = function() return(private$r_gammas),
        invariant = function() return(private$r_invariant),
        model = function() return(private$r_model)



    ),

    private = list(

        r_Q = NULL,
        r_pi_tcag = NULL,
        r_U = NULL,
        r_Ui = NULL,
        r_L = NULL,
        r_gammas = NULL,
        r_invariant = NULL,
        r_model = NULL

    ),

    lock_class = TRUE

)









# .... partial functions -------


#' Check sub_* function arguments for validity.
#'
#' @noRd
#'
sub_arg_checks <- function(mu, mod_name,
                           pi_tcag = NULL, alpha_1 = NULL, alpha_2 = NULL,
                           beta = NULL, gamma_shape = NULL, gamma_k = NULL,
                           invariant = NULL, lambda = NULL, alpha = NULL,
                           kappa = NULL, abcdef = NULL, Q = NULL) {

    # Vector parameters:
    if (!is.null(pi_tcag) && !(is_type(pi_tcag, c("integer", "numeric"), 4) &&
                               all(pi_tcag >= 0) && any(pi_tcag > 0))) {
        err_msg(paste0("sub_", mod_name), "pi_tcag",
                "a length-4 numeric vector where at least one number is > 0 and all",
                "are >= 0")
    }
    if (!is.null(abcdef) && !(is_type(abcdef, c("integer", "numeric"), 6) &&
                              all(abcdef >= 0))) {
        err_msg(paste0("sub_", mod_name), "abcdef",
                "a length-6 numeric vector where all numbers are >= 0")
    }

    # UNREST matrix:
    if (!is.null(Q) && !(is_type(Q, "matrix") && nrow(Q) == 4 && ncol(Q) == 4 &&
                         is.numeric(Q) &&
                         all(c(Q[lower.tri(Q)], Q[upper.tri(Q)]) >= 0))) {
        err_msg(paste0("sub_", mod_name), "Q",
                "a 4x4 numeric matrix where all non-diagonal elements are >= 0")
    }

    # Single-number parameters:
    if (!is.null(mu) && !(single_number(mu) && mu > 0)) {
        err_msg(paste0("sub_", mod_name), "mu", "NULL or a single number > 0")
    }
    if (!is.null(alpha_1) && !single_number(alpha_1, 0)) {
        err_msg(paste0("sub_", mod_name), "alpha_1", "a single number >= 0")
    }
    if (!is.null(alpha_2) && !single_number(alpha_2, 0)) {
        err_msg(paste0("sub_", mod_name), "alpha_2", "a single number >= 0")
    }
    if (!is.null(beta) && !single_number(beta, 0)) {
        err_msg(paste0("sub_", mod_name), "beta", "a single number >= 0")
    }
    if (!is.null(lambda) && !single_number(lambda, 0)) {
        err_msg(paste0("sub_", mod_name), "lambda", "a single number >= 0")
    }
    if (!is.null(alpha) && !single_number(alpha, 0)) {
        err_msg(paste0("sub_", mod_name), "alpha", "a single number >= 0")
    }
    if (!is.null(kappa) && !single_number(kappa, 0)) {
        err_msg(paste0("sub_", mod_name), "kappa", "a single number >= 0")
    }


    # Site-heterogeneity parameters:
    if (!is.null(gamma_shape) && !(single_number(gamma_shape) && gamma_shape > 0)) {
        err_msg(paste0("sub_", mod_name), "gamma_shape",
                "NULL or a single number > 0")
    }
    if (!single_integer(gamma_k, 2, 255)) {
        err_msg(paste0("sub_", mod_name), "gamma_k",
                "a single integer in range [2, 255]")
    }
    if (!single_number(invariant, 0) || invariant >= 1) {
        err_msg(paste0("sub_", mod_name), "invariant",
                "a single number >= 0 and < 1")
    }


    invisible(NULL)


}


# ... sub_models docs -----
#' Construct necessary information for substitution models.
#'
#' For a more detailed explanation, see `vignette("sub-models")`.
#'
#'
#' @name sub_models
#'
#' @seealso \code{\link{create_haplotypes}}
#'
#' @return A `sub_info` object, which is an R6 class that wraps the info needed for
#' the `create_haplotypes` function.
#' It does not allow the user to directly manipulate the info inside, as that
#' should be done using the `sub_models` functions.
#' You can use the following methods from the class to view information:
#' \describe{
#'     \item{`Q()`}{View a list of substitution rate matrices,
#'         one for each Gamma category.}
#'     \item{`pi_tcag()`}{View the equilibrium nucleotide frequencies.}
#'     \item{`gammas()`}{View the discrete Gamma-class values.}
#'     \item{`invariant()`}{View the proportion of invariant sites.}
#'     \item{`model()`}{View the substitution model.}
#'     \item{`U()`}{View list of the `U` matrices (one matrix per Gamma category)
#'         used for calculating transition-probability matrices.
#'         This is empty for UNREST models.}
#'     \item{`Ui()`}{View list of the `U^-1` matrices (one matrix per Gamma category)
#'         used for calculating transition-probability matrices.
#'         This is empty for UNREST models.}
#'     \item{`L()`}{View list of the lambda vectors (one vector per Gamma category)
#'         used for calculating transition-probability matrices.
#'         This is empty for UNREST models.}
#' }
#'
#'
#'
#' @examples
#' # Same substitution rate for all types:
#' obj_JC69 <- sub_JC69(lambda = 0.1)
#'
#' # Transitions 2x more likely than transversions:
#' obj_K80 <- sub_K80(alpha = 0.2, beta = 0.1)
#'
#' # Incorporating equilibrium frequencies:
#' obj_HKY85 <- sub_HKY85(pi_tcag = c(0.1, 0.2, 0.3, 0.4),
#'                        alpha = 0.2, beta = 0.1)
#'
#' # 10-category Gamma distribution for among-site variability:
#' obj_K80 <- sub_K80(alpha = 0.2, beta = 0.1,
#'                    gamma_shape = 1, gamma_k = 10)
#'
#' # Invariant sites:
#' obj_K80 <- sub_K80(alpha = 0.2, beta = 0.1,
#'                    invariant = 0.25)
#'
NULL




#' @describeIn sub_models JC69 model.
#'
#' @param lambda Substitution rate for all possible substitutions.
#'
#' @export
#'
#'
sub_JC69 <- function(lambda, mu = 1, gamma_shape = NULL, gamma_k = 5, invariant = 0) {

    sub_arg_checks(mu, "JC69", lambda = lambda,
                   gamma_shape = gamma_shape, gamma_k = gamma_k, invariant = invariant)

    pi_tcag <- rep(0.25, 4)
    lambda <- lambda * 4;  # bc it's being multiplied by pi_tcag

    out <- sub_TN93__(mu, pi_tcag, lambda, lambda, lambda,
                    gamma_shape, gamma_k, invariant, "JC69")

    return(out)

}

#' @describeIn sub_models K80 model.
#'
#' @param alpha Substitution rate for transitions.
#'
#' @export
#'
sub_K80 <- function(alpha, beta, mu = 1, gamma_shape = NULL, gamma_k = 5, invariant = 0) {

    sub_arg_checks(mu, "K80", alpha = alpha, beta = beta,
                   gamma_shape = gamma_shape, gamma_k = gamma_k, invariant = invariant)

    pi_tcag <- rep(0.25, 4)
    alpha <- alpha * 4;  # bc they're being multiplied by pi_tcag
    beta <- beta * 4;  # bc they're being multiplied by pi_tcag

    out <- sub_TN93__(mu, pi_tcag, alpha, alpha, beta,
                    gamma_shape, gamma_k, invariant, "K80")

    return(out)

}

#' @describeIn sub_models F81 model.
#'
#'
#' @export
#'
sub_F81 <- function(pi_tcag, mu = 1, gamma_shape = NULL, gamma_k = 5, invariant = 0) {

    sub_arg_checks(mu, "F81", pi_tcag = pi_tcag,
                   gamma_shape = gamma_shape, gamma_k = gamma_k, invariant = invariant)

    out <- sub_TN93__(mu, pi_tcag, 1, 1, 1,
                    gamma_shape, gamma_k, invariant, "F81")

    return(out)

}

#' @describeIn sub_models HKY85 model.
#'
#'
#'
#' @export
#'
sub_HKY85 <- function(pi_tcag, alpha, beta,
                      mu = 1, gamma_shape = NULL, gamma_k = 5, invariant = 0) {

    sub_arg_checks(mu, "HKY85", pi_tcag = pi_tcag, alpha = alpha, beta = beta,
                   gamma_shape = gamma_shape, gamma_k = gamma_k, invariant = invariant)

    out <- sub_TN93__(mu, pi_tcag, alpha, alpha, beta,
                    gamma_shape, gamma_k, invariant, "HKY85")

    return(out)

}

#' @describeIn sub_models F84 model.
#'
#'
#' @param kappa The transition/transversion rate ratio.
#'
#' @export
#'
sub_F84 <- function(pi_tcag, beta, kappa,
                    mu = 1, gamma_shape = NULL, gamma_k = 5, invariant = 0) {


    sub_arg_checks(mu, "F84", pi_tcag = pi_tcag, beta = beta, kappa = kappa,
                   gamma_shape = gamma_shape, gamma_k = gamma_k, invariant = invariant)

    pi_y = pi_tcag[1] + pi_tcag[2]
    pi_r = pi_tcag[3] + pi_tcag[4]

    alpha_1 = (1 + kappa / pi_y) * beta
    alpha_2 = (1 + kappa / pi_r) * beta

    out <- sub_TN93__(mu, pi_tcag, alpha_1, alpha_2, beta,
                    gamma_shape, gamma_k, invariant, "F84")

    return(out)

}



#' @describeIn sub_models TN93 model.
#'
#' @param pi_tcag Vector of length 4 indicating the equilibrium distributions of
#'     T, C, A, and G respectively. Values must be >= 0, and
#'     they are forced to sum to 1.
#' @param alpha_1 Substitution rate for T <-> C transition.
#' @param alpha_2 Substitution rate for A <-> G transition.
#' @param beta Substitution rate for transversions.
#' @param mu Total rate of substitutions. Defaults to `1`, which makes branch lengths
#'     in units of substitutions per site. Passing `NULL` results in no scaling.
#' @param gamma_shape Numeric shape parameter for discrete Gamma distribution used for
#'     among-site variability. Values must be greater than zero.
#'     If this parameter is `NULL`, among-site variability is not included.
#'     Defaults to `NULL`.
#' @param gamma_k The number of categories to split the discrete Gamma distribution
#'     into. Values must be an integer in the range `[2,255]`.
#'     This argument is ignored if `gamma_shape` is `NA`.
#'     Defaults to `5`.
#' @param invariant Proportion of sites that are invariant.
#'     Values must be in the range `[0,1)`.
#'     Defaults to `0`.
#'
#' @export
#'
sub_TN93 <- function(pi_tcag, alpha_1, alpha_2, beta,
                     mu = 1, gamma_shape = NULL, gamma_k = 5, invariant = 0) {

    out <- sub_TN93__(mu, pi_tcag, alpha_1, alpha_2, beta, gamma_shape, gamma_k,
                            invariant, "TN93")

    return(out)

}



# .... full functions -----
# (The above functions are simply special cases of `sub_TN93` and use that
#  function internally.)



#' Inner function that does most of the work for TN93 and its special cases
#'
#' @noRd
#'
sub_TN93__ <- function(mu, pi_tcag, alpha_1, alpha_2, beta,
                       gamma_shape, gamma_k, invariant, model) {

    sub_arg_checks(mu, "TN93", pi_tcag = pi_tcag,
                   alpha_1 = alpha_1, alpha_2 = alpha_2, beta = beta,
                   gamma_shape = gamma_shape, gamma_k = gamma_k, invariant = invariant)
    if (is.null(gamma_shape)) gamma_shape <- 0

    if (!is_type(model, "character", 1L)) {
        stop("\nINNER ERROR: arg `model` to sub_TN93__ is not a single string.")
    }
    if (is.null(mu)) mu <- -1

    info_list <- sub_TN93_cpp(mu, pi_tcag, alpha_1, alpha_2, beta, gamma_shape, gamma_k,
                              invariant)
    info_list[["model"]] <- model

    out <- sub_info$new(info_list)

    return(out)

}



#' @describeIn sub_models GTR model.
#'
#' @param abcdef A vector of length 6 that contains the off-diagonal elements
#'     for the substitution rate matrix.
#'     See `vignette("sub-models")` for how the values are ordered in the matrix.
#'
#' @export
#'
sub_GTR <- function(pi_tcag, abcdef, mu = 1, gamma_shape = NULL, gamma_k = 5, invariant = 0) {

    sub_arg_checks(mu, "GTR", pi_tcag = pi_tcag, abcdef = abcdef,
                   gamma_shape = gamma_shape, gamma_k = gamma_k, invariant = invariant)
    if (is.null(gamma_shape)) gamma_shape <- 0

    if (is.null(mu)) mu <- -1

    info_list <- sub_GTR_cpp(mu, pi_tcag, abcdef, gamma_shape, gamma_k, invariant)

    out <- sub_info$new(info_list)

    return(out)

}

#' @describeIn sub_models UNREST model.
#'
#'
#' @param Q Matrix of substitution rates for "T", "C", "A", and "G", respectively.
#'     Item `Q[i,j]` is the rate of substitution from nucleotide `i` to nucleotide `j`.
#'     Do not include indel rates here!
#'     Values on the diagonal are calculated inside the function so are ignored.
#'
#' @export
#'
#'
sub_UNREST <- function(Q, mu = 1, gamma_shape = NULL, gamma_k = 5, invariant = 0) {

    sub_arg_checks(mu, "UNREST", Q = Q,
                   gamma_shape = gamma_shape, gamma_k = gamma_k, invariant = invariant)
    if (is.null(gamma_shape)) gamma_shape <- 0

    if (is.null(mu)) mu <- -1

    info_list <- sub_UNREST_cpp(mu, Q, gamma_shape, gamma_k, invariant)

    out <- sub_info$new(info_list)

    return(out)

}






# INDELS -----


# ... indel_info class -------



#' An R6 class representing information for indel rates.
#'
#'
#' This class should NEVER be created using `indel_info$new`.
#' Only use function `indels`.
#' That's why I'm not exporting it.
#'
#' @noRd
#'
#'
#' @importFrom R6 R6Class
#'
indel_info <- R6Class(

    "indel_info",

    public = list(

        initialize = function(rates) {

            extra_msg <- paste(" Please only create these objects using the indels",
                               "function, NOT using indel_info$new().")
            if (!inherits(rates, "numeric") || any(rates < 0)) {
                stop("\nWhen initializing an indel_info object, you need to use ",
                     "a numeric vector of values >= 0.", extra_msg, call. = FALSE)
            }

            private$r_rates <- rates
        },

        print = function(...) {

            digits <- max(3, getOption("digits") - 3)

            cat("< Indel rates >\n")
            cat(sprintf(sprintf("# Total rate = %%.%ig\n", digits), sum(private$r_rates)))
            cat(sprintf("# Max length = %i\n", length(private$r_rates)))

            invisible(self)

        },

        rates = function() return(private$r_rates)

    ),

    private = list(

        r_rates = NULL

    ),

    lock_class = TRUE

)





# ... main function -----


#' Insertions and deletions (indels) specification
#'
#' Construct necessary information for insertions and deletions (indels) that will
#' be used in `create_haplotypes`.
#'
#' All indels require the `rate` parameter, which specifies
#' the overall indels rate among all lengths.
#' The `rate` parameter is ultimately combined with a vector of relative rates among
#' the different lengths of indels from 1 to the maximum possible length.
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
#'         `{L * max_length / (max_length - L + 1)}^(-a)`.
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
#'
#' @param rate Single number specifying the overall indel rate among all lengths.
#' @param max_length Maximum length of indels. Defaults to `10`.
#' @param a Extra parameter necessary for generating rates from a Lavalette distribution.
#'     See Details for more info. Defaults to `NULL`.
#' @param rel_rates A numeric vector of relative rates for each indel length
#'     from 1 to the maximum length.
#'     If provided, all arguments other than `rate` are ignored.
#'     Defaults to `NULL`.
#'
#'
#' @references
#' Albers, C. A., G. Lunter, D. G. MacArthur, G. McVean, W. H. Ouwehand, and R. Durbin.
#' 2011. Dindel: accurate indel calls from short-read data. Genome Research 21:961–973.
#'
#' Fletcher, W., and Z. Yang. 2009. INDELible: a flexible simulator of
#' biological sequence evolution. Molecular Biology and Evolution 26:1879–1888.
#'
#' @export
#'
#' @return An `indel_info` object, which is an R6 class that wraps the info needed for
#' the `create_haplotypes` function.
#' It does not allow the user to directly manipulate the info inside, as that
#' should be done using this function.
#' You can use the `rates()` method to view the indel rates by size.
#'
#' @examples
#' # relative rates are proportional to `exp(-L)` for indel
#' # length `L` from 1 to 5:
#' indel_rates1 <- indels(0.1, max_length = 5)
#'
#' # relative rates are proportional to Lavalette distribution
#' # for length from 1 to 10:
#' indel_rates2 <- indels(0.2, max_length = 10, a = 1.1)
#'
#' # relative rates are all the same for lengths from 1 to 100:
#' indel_rates3 <- indels(0.2, rel_rates = rep(1, 100))
#'
#'
#'
indels <- function(rate,
                   max_length = 10,
                   a = NULL,
                   rel_rates = NULL) {

    # Checking types:
    if (!single_number(rate) || rate <= 0) {
        err_msg("indels", "rate", "a single number > 0")
    }
    if (!single_integer(max_length, 1, 1e6)) {
        err_msg("indels", "max_length", "a single integer in range [1, 1e6]")
    }
    if (!is.null(a) && !single_number(a, 0)) {
        err_msg("indels", "a", "NULL or a single number >= 0")
    }
    if (!is.null(rel_rates) && !positive_vector(rel_rates)) {
        err_msg("indels", "rel_rates", "NULL or a vector of positive numbers that ",
                "sums to > 0")
    }

    # Generate relative rates if not directly specified:
    if (is.null(rel_rates)) {
        if (!is.null(a)) {
            L <- 1:(max_length)
            rel_rates <- {(L * max_length) / (max_length - L + 1)}^(-a)
        } else {
            rel_rates <- exp(-1 * 1:(max_length))
        }
    }

    # So relative rates sum to 1:
    rel_rates <- rel_rates / sum(rel_rates)

    # Absolute rates:
    rates <- rel_rates * rate

    rates_obj <- indel_info$new(rates)

    return(rates_obj)
}

