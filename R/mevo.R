#' Check validity of input parameters for substitions and indels.
#'
#'
#' @noRd
#'
validate_sub_pars <- function(par_list) {

    # Defaults err_msg fxn doesn't quite work
    err_msg_ <- "\nThe parameter `%s` used for substititions should be %s."

    for (x in c("lambda", "alpha", "beta", "alpha_1", "alpha_2", "kappa")) {
        z <- par_list[[x]]
        if (!is.null(z) && !single_number(z, 0)) {
            stop(sprintf(err_msg_, x, "a single number >= 0"), call. = FALSE)
        }
    }
    if (!is.null(par_list[["pi_tcag"]]) &&
        (!is_type(par_list[["pi_tcag"]], "numeric", 4) ||
         any(par_list[["pi_tcag"]] < 0) ||
         all(par_list[["pi_tcag"]] == 0))) {
        stop(sprintf(err_msg_, "pi_tcag",
                     paste("numeric vector of length 4, where",
                           "no values should be < 0 and",
                           "at least one value should be > 0")), call. = FALSE)
    }
    if (!is.null(par_list[["abcdef"]]) &&
        (!is_type(par_list[["abcdef"]], "numeric", 6) ||
         any(par_list[["abcdef"]] < 0) ||
         all(par_list[["abcdef"]] == 0))) {
        stop(sprintf(err_msg_, "abcdef",
                     paste("numeric vector of length 6, where",
                           "no values should be < 0 and",
                           "at least one value should be > 0")), call. = FALSE)
    }
    if (!is.null(par_list[["Q"]]) && (!is_type(par_list[["Q"]], "matrix") ||
                        !identical(dim(par_list[["Q"]]), c(4L, 4L)) ||
                        any(par_list[["Q"]] < 0))) {
        stop(sprintf(err_msg_, "Q", "a 4x4 matrix, where no values are < 0"),
             call. = FALSE)
    }

    invisible(NULL)
}


#' Make substitution rate matrix (Q) and vector of equilibrium frequencies (pis).
#'
#' @inheritParams make_mevo
#'
#' @noRd
#'
#'
substitutions <- function(sub) {

    err_msg_ <- paste("\nNot all required parameters provided for the given",
                     "substitution model (\"%s\").",
                     "See `vignette(\"sub-model\")` for each model's required",
                     "parameter(s).")

    if (sub$model == "TN93") {
        if (any(! c("pi_tcag", "alpha_1", "alpha_2", "beta") %in% names(sub))) {
            stop(sprintf(err_msg_, sub$model), call. = FALSE)
        }
        if (any(sub$pi_tcag < 0) || sum(sub$pi_tcag <= 0)) {
            stop("\nNo values of `pi_tcag` should be < 0, and the vector needs at ",
                 "least one value > 0.", call. = FALSE)
        }
        Q <- TN93_rate_matrix(sub$pi_tcag, sub$alpha_1, sub$alpha_2,
                              sub$beta)
        pi_tcag <- sub$pi_tcag
    } else if (sub$model == "JC69") {
        if (any(! c("lambda") %in% names(sub))) {
            stop(sprintf(err_msg_, sub$model), call. = FALSE)
        }
        Q <- JC69_rate_matrix(sub$lambda)
        pi_tcag <- rep(0.25, 4)
    } else if (sub$model == "K80") {
        if (any(! c("alpha", "beta") %in% names(sub))) {
            stop(sprintf(err_msg_, sub$model), call. = FALSE)
        }
        Q <- K80_rate_matrix(sub$alpha, sub$beta);
        pi_tcag <- rep(0.25, 4)
    } else if (sub$model == "F81") {
        if (any(! c("pi_tcag") %in% names(sub))) {
            stop(sprintf(err_msg_, sub$model), call. = FALSE)
        }
        Q <- F81_rate_matrix(sub$pi_tcag);
        pi_tcag <- sub$pi_tcag
    } else if (sub$model == "HKY85") {
        if (any(! c("pi_tcag", "alpha", "beta") %in% names(sub))) {
            stop(sprintf(err_msg_, sub$model), call. = FALSE)
        }
        Q <- HKY85_rate_matrix(sub$pi_tcag, sub$alpha, sub$beta)
        pi_tcag <- sub$pi_tcag
    } else if (sub$model == "F84") {
        if (any(! c("pi_tcag", "beta", "kappa") %in% names(sub))) {
            stop(sprintf(err_msg_, sub$model), call. = FALSE)
        }
        Q <- F84_rate_matrix(sub$pi_tcag, sub$beta, sub$kappa)
        pi_tcag <- sub$pi_tcag
    } else if (sub$model == "GTR") {
        if (any(! c("pi_tcag", "abcdef") %in% names(sub))) {
            stop(sprintf(err_msg_, sub$model), call. = FALSE)
        }
        Q <- GTR_rate_matrix(sub$pi_tcag, sub$abcdef)
        pi_tcag <- sub$pi_tcag
    } else if (sub$model == "UNREST") {
        if (any(! c("Q") %in% names(sub))) {
            stop(sprintf(err_msg_, sub$model), call. = FALSE)
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
        err_msg_ <- paste0("\nWhen specifying `", which_type, "` in `make_mevo`, ")
        if (!single_number(indel$rate, 0)) {
            err_msg("make_mevo", which_type, "a list with a \"rate\" field,",
                    "and that field must be a single number >= 0")
        }
        if (indel$rate == 0) return(numeric(0))
        names_ <- sort(names(indel))
        if (identical(names_, sort(c("rate", "max_length")))) {
            if (!single_integer(indel$max_length, 1)) {
                stop(err_msg_, "the `max_length` field within, if provided, must ",
                     "be a single whole number >= 1.", call. = FALSE)
            }
            rel_rates <- exp(-1 * 1:(indel$max_length))
        } else if (identical(names_, sort(c("rate", "max_length", "a")))) {
            if (!single_integer(indel$max_length, 1)) {
                stop(err_msg_, "the `max_length` field within, if provided, must ",
                     "be a single whole number >= 1.", call. = FALSE)
            }
            if (!single_number(indel$a, 0)) {
                stop(err_msg_, "the `a` field within, if provided, must ",
                     "be a single number >= 0.", call. = FALSE)
            }
            L <- 1:(indel$max_length)
            rel_rates <- {(L * indel$max_length) / (indel$max_length - L + 1)}^(-indel$a)
        } else if (identical(names_, sort(c("rate", "rel_rates")))) {
            if (!is_type(indel$rel_rates, "numeric") ||
                 any(indel$rel_rates < 0) ||
                 all(indel$rel_rates == 0)) {
                stop(sprintf(err_msg_, "pi_tcag",
                             paste("numeric vector of length 4, where",
                                   "no values should be < 0 and",
                                   "at least one value should be > 0")), call. = FALSE)
                stop(err_msg_, "the `rel_rates` field within, if provided, must ",
                     "be a numeric vector, where no values should be < 0 and ",
                     "at least one value should be > 0.", call. = FALSE)
            }
            rel_rates <- indel$rel_rates
        } else {
            err_msg("make_mevo", which_type, "a list with names that coincide",
                    "with one of the methods in the \"Indels\" section within",
                    "`?make_mevo`. Note that extra names return an error.")
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

            err_msg_ <- paste("\nThe `mats` field inside the `site_var`",
                             "argument to the `make_mevo` function needs to",
                             "be a list of matrices.")
            if (!inherits(site_var$mats, "list")) {
                stop(err_msg_, call. = FALSE)
            } else if (!all(sapply(site_var$mats, inherits,
                                   what = "matrix"))) {
                stop(err_msg_, call. = FALSE)
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
    validate_sub_pars(sub)
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







#' Process one segregating-sites matrix from a coalescent simulator with ms-style output.
#'
#' @param mat The matrix to process.
#' @param seq_size The number of bp in the sequence associated with the input string.
#'
#' @noRd
#'
process_coal_obj_sites <- function(mat, seq_size) {

    err_msg_ <- paste("\nPositions in one or more segregating-sites matrices %s.",
                      "They are derived from column names, so check those for nonsense.")

    # Dealing with coala objects:
    if (!is.null(mat$snps)) mat <- mat$snps

    pos <- as.numeric(colnames(mat))

    if (any(is.na(pos))) {
        stop(sprintf(err_msg_, "are producing NAs when trying to coerce them to numbers"),
             call. = FALSE)
    } else if (length(pos) != ncol(mat)) {
        stop(sprintf(err_msg_, paste("are not the same length as the number of rows",
                                     "in the matrix")), call. = FALSE)
    }

    new_mat <- cbind(pos, t(mat))
    rownames(new_mat) <- colnames(new_mat) <- NULL

    return(new_mat)
}



#' Read segregating-site info from a coalescent object from scrm or coala.
#'
#' @inheritParams make_phylo_info
#'
#' @return An XPtr to the info needed from the phylogenies to do the sequence
#'     simulations.
#'
#' @noRd
#'
#'
coal_obj_sites <- function(coal_obj) {

    # Check for coal_obj being a list and having a `seg_sites` field
    if (!inherits(coal_obj, "list") || is.null(coal_obj$seg_sites)) {
        err_msg("create_variants", "method_info",
                "a list with a `seg_sites` field present (when method = \"coal_sites\").",
                "This must also be true of the inner list if also providing a",
                "`names` field.")
    }

    sites <- coal_obj$seg_sites

    sites_mats <- lapply(sites, process_coal_obj_sites)

    return(sites_mats)
}


#' Check validity of position columns in segregating-sites matrices.
#'
#' @noRd
#'
fill_coal_mat_pos <- function(sites_mats, seq_sizes) {

    if (length(sites_mats) != length(seq_sizes)) {
        err_msg("create_variants", "method_info",
                "a list with `seg_sites` field that contains the same number of",
                "matrices as the number of sequences.")
    }

    for (i in 1:length(sites_mats)) {
        if (nrow(sites_mats[[i]]) == 0) next;
        if (all(sites_mats[[i]][,1] < 1 && sites_mats[[i]][,1] > 0)) {
            # Converting to integer positions (0-based):
            sites_mats[[i]][,1]  <- as.integer(sites_mats[[i]][,1] * seq_sizes[i]);
        } else if (all(sites_mats[[i]][,1] < seq_sizes[i] && sites_mats[[i]][,1] >= 0)) {
            # Keeping them in 0-based indices:
            sites_mats[[i]][,1] = as.integer(sites_mats[[i]][,1]);
        } else if (all(sites_mats[[i]][,1] <= seq_sizes[i] && sites_mats[[i]][,1] >= 1)) {
            # Converting to 0-based indices:
            sites_mats[[i]][,1] = as.integer(sites_mats[[i]][,1]) - 1;
        } else {
            stop("\nYour `method_info` argument to `create_variants` is causing problems: ",
                 "Positions in one or more segregating-sites matrices ",
                 "are not obviously from either a finite- or infinite-sites model. ",
                 "The former should have positions in the range ",
                 "[0, sequence length - 1] or [1, sequence length], ",
                 "the latter in (0,1).", call. = FALSE)
        }

    }

    return(sites_mats)

}


#' Create variants from segregating-site info from coalescent simulations.
#'
#'
#' @noRd
#'
read_coal_sites <- function(method_info, reference, mevo_obj, n_cores, show_progress) {

    if (inherits(method_info, "list") && !is.null(method_info$names) &&
        !is.null(method_info$info)) {
        var_names <- method_info$names
        method_info <- method_info$info
    } else {
        var_names <- character(0)
    }

    if (inherits(method_info, "list")) {
        sites_mats <- coal_obj_sites(method_info)
    } else if (is_type(method_info, "character", 1)) {
        sites_mats <- coal_file_sites(method_info)
        # Revert back to list (from arma::field which adds dims)
        dim(sites_mats) <- NULL
    } else {
        err_msg("create_variants", "method_info", "a single string or a list",
                "(for method = \"coal_sites\")")
    }


    if (length(var_names) == 0) {
        var_names <- sprintf("var%i", 0:(ncol(sites_mats[[1]])-2))
    }
    if (length(unique(sapply(sites_mats, ncol))) != 1) {
        err_msg("create_variants", "method_info", "a list containing matrices with",
                "the same number of rows, if `method` is \"coal_sites\"")
    }
    if (unique(sapply(sites_mats, ncol))[1] != (length(var_names) + 1)) {
        err_msg("create_variants", "method_info",
                "(when `method` is \"coal_sites\" and when providing a names field in",
                "`method_info`)",
                "a list containing a `names` field",
                "with the same number of items as the number of rows in each",
                "segregating-sites matrix")
    }

    seq_sizes <- reference$sizes()

    # Fill and check the position column in sites_mats
    sites_mats <- fill_coal_mat_pos(sites_mats, seq_sizes)

    variants_ptr <- add_coal_sites_cpp(reference$genome,
                                       var_names,
                                       sites_mats,
                                       mevo_obj$Q,
                                       mevo_obj$pi_tcag,
                                       mevo_obj$insertion_rates,
                                       mevo_obj$deletion_rates,
                                       n_cores,
                                       show_progress)

    return(variants_ptr)

}


