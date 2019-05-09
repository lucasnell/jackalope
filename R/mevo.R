
#' ========================================================================`
#'
#' This file stores functions for molecular-evolution info.
#' This includes those for substitutions, indels,
#' among-site variability in mutation rates, and
#' creating `mevo` objects that organize molecular-evolution info.
#'
#' ========================================================================`


# substitutions -----


#' Print method for `sub_model_info` objects, output from `sub_models` functions.
#'
#' I added this mostly to make users less likely to edit it manually and to
#' give context to the output.
#'
#' @noRd
#' @export
#'
print.sub_model_info <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("< Substitution information >\n")
    fmt <- paste0("%.", digits, "f")

    cat("# Equilibrium densities:\n")
    cat("  ", sprintf(fmt, x$pi_tcag), "\n")

    cat("# Substitution rate matrix:\n")
    prmatrix(x$Q, digits = digits,
             rowlab = paste("  ", c("T", "C", "A", "G")),
             collab = c("T", "C", "A", "G"))

    cat("(View rate matrix in the \"Q\" field)\n")
    cat("(View equil. densities in the \"pi_tcag\" field)")

    invisible(NULL)
}





# indels -----


#' Insertions and deletions (indels) specification
#'
#' Construct necessary information for insertions and deletions (indels) that will
#' be used in `create_variants`.
#'
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
#' @return An `indel_rates` object, which is just a wrapper around a numeric vector.
#' You can access the rates vector for `indel_rates` object `x` by running
#' `as.numeric(x)`.
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

    class(rates) <- "indel_rates"

    return(rates)
}

#' Print method for indel_rates objects.
#'
#' I added this mostly to make sure a giant vector doesn't ever print.
#'
#' @noRd
#' @export
#'
print.indel_rates <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("< Indel rates vector >\n")
    cat(sprintf("  * Total rate = %.3g\n", sum(x)))
    cat(sprintf("  * Max length = %i\n", length(x)))
    cat("(View raw data using `as.numeric`)\n")
    invisible(NULL)
}



# among-site var. -----


#' Specify variation in mutation rates among sites
#'
#' Construct necessary information for among-site variation in mutation rates that will
#' be used in `create_variants`.
#'
#'
#' A site's deviance from the average mutation rate is determined by its
#' "gamma distance".
#' A site's overall mutation rate is the mutation rate for that nucleotide
#' (substitution + indel) multiplied by the site's gamma distance.
#' There are two options for specifying gamma distances:
#' \enumerate{
#'     \item Generate gamma distances from a Gamma distribution.
#'         This method will be used if the `shape` and `region_size` arguments
#'         are provided.
#'         If the `mats` argument is also provided, this method will NOT be used.
#'         See argument descriptions for more info.
#'     \item Manually input matrices that specify the gamma distance and end points
#'         for regions each gamma distance refers to.
#'         This method will be used if the `mats` argument is provided.
#'         See argument descriptions for more info.
#' }
#'
#' @param reference A \code{ref_genome} object from which you will eventually
#'     generate variants.
#' @param shape Shape parameter for the Gamma distribution that generates gamma distances,
#'     The variance of the distribution is `1 / shape`, and its mean is fixed to 1.
#'     Values `<= 0` are not allowed.
#'     Defaults to `NULL`.
#' @param region_size Size of regions to break the genome into,
#'     where all sites within a region have the same gamma distance.
#'     Defaults to `NULL`.
#' @param invariant Proportion of regions that are invariant.
#'     Must be in range `[0,1)`.
#'     Defaults to `0`.
#' @param mats List of matrices, one for each sequence in the genome.
#'     Each matrix should have two columns.
#'     The first should contain the end points for each region.
#'     The second should contain the gamma distances for each region.
#'     Note that if gamma distances don't have a mean (weighted by
#'     sequence length for each gamma-distance value) equal to 1,
#'     you're essentially changing the overall mutation rate.
#'     If this argument is provided, `shape` and `region_size` are ignored.
#'     Defaults to `NULL`.
#' @param out_prefix String specifying the file name prefix for an output BED file that
#'     will be generated by this function and that will specify the
#'     gamma distances for each region.
#'     If `NULL`, no output file is produced.
#'     Defaults to `NULL`.
#' @inheritParams write_fasta
#'
#'
#' @export
#'
#' @return A `site_var_mats` object, which is a wrapper around a list of matrices,
#' one for each sequence in the reference genome.
#' Although the print method is different, you can otherwise treat these objects
#' the same as you would a list (e.g., `x[[1]]`, `x[1:2]`, `length(x)`).
#'
#'
#' @examples
#' ref <- create_genome(3, 100)
#' # generating from Gamma distribution
#' gamma_mats <- site_var(ref, shape = 0.5,
#'                        region_size = 5)
#' # with custom matrices
#' gamma_mats <- site_var(ref,
#'                        mats = replicate(3,
#'                            cbind(seq(10, 100, 10),
#'                            rgamma(10, 0.9)),
#'                            simplify = FALSE))
#'
site_var <- function(reference,
                     shape = NULL,
                     region_size = NULL,
                     invariant = 0,
                     mats = NULL,
                     out_prefix = NULL,
                     compress = FALSE,
                     comp_method = "bgzip") {

    # ---------*
    # Checking types:
    # ---------*
    if (!inherits(reference, "ref_genome")) {
        err_msg("site_var", "reference", "a \"ref_genome\" object")
    }
    if (!is.null(shape) && (!single_number(shape) || shape <= 0)) {
        err_msg("site_var", "shape", "NULL or a single number > 0")
    }
    if (!is.null(region_size) && !single_integer(region_size, 1)) {
        err_msg("site_var", "region_size", "NULL or a single integer >= 1")
    }
    if (!is.null(invariant) && (!single_number(invariant, 0) || invariant >= 1)) {
        err_msg("site_var", "invariant", "NULL or a single number in range [0,1)")
    }
    if (!is.null(mats) && (!is_type(mats, "list") ||
                           !all(sapply(mats, inherits, what = "matrix")))) {
        err_msg("site_var", "mats", "NULL or a list of matrices")
    }
    if (!is.null(out_prefix) && !is_type(out_prefix, "character", 1)) {
        err_msg("site_var", "out_prefix", "NULL or a single string")
    }
    if (!is_type(compress, "logical", 1) && !single_integer(compress, 1, 9)) {
        err_msg("site_var", "compress", "a single logical or integer from 1 to 9")
    }
    if (is_type(compress, "logical", 1) && compress) compress <- 6 # default compression
    if (is_type(compress, "logical", 1) && !compress) compress <- 0 # no compression
    if (!is_type(comp_method, "character", 1) || !comp_method %in% c("gzip", "bgzip")) {
        err_msg("site_var", "comp_method", "\"gzip\" or \"bgzip\"")
    }

    # Checking for other nonsense:
    if ((is.null(shape) && !is.null(region_size)) ||
        (!is.null(shape) && is.null(region_size))) {
        stop("\nIn the `site_var` function, if you provide an input to the ",
             "`region_size` argument, you must also provide one to the `shape` ",
             "argument, and vice versa.", call. = FALSE)
    }

    # ---------*
    # Making the matrices:
    # ---------*
    seq_sizes <- reference$sizes()
    if (!is.null(mats)) {
        if (length(mats) != length(seq_sizes)) {
            err_msg("site_var", "mats", "NULL or a list of matrices the same length",
                    "as the number of sequences in the reference genome.")
        }
    } else {
        if (is.null(shape) || is.null(region_size)) {
            stop("\nIn function `site_var`, if you don't provide a `mats` argument, ",
                 "you need to provide both the `shape` and `region_size` arguments.",
                 call. = FALSE)
        }
        mats <- make_gamma_mats(seq_sizes,
                                region_size_ = region_size,
                                shape = shape,
                                invariant = invariant)
        dim(mats) <- NULL # so it's just a list now
    }

    # Check matrices for proper end points and # columns:
    check_gamma_mats(mats, seq_sizes)

    # ---------*
    # Writing to BED file if desired:
    # ---------*
    if (!is.null(out_prefix)) {
        seq_names <- reference$names()
        write_bed(out_prefix, mats, seq_names, compress, comp_method)
    }

    class(mats) <- "site_var_mats"

    return(mats)
}


#' Print output from `site_var`
#'
#' I added this mostly to make sure a giant list doesn't print.
#'
#' @noRd
#' @export
#'
print.site_var_mats <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("< Site variability matrices >\n")
    cat(sprintf("  * Number of sequences = %i\n", length(x)))
    cat(sprintf("  * Number of regions = %i\n", sum(sapply(x, nrow))))
    cat("(View raw data the same as you would a list)\n")
    invisible(NULL)
}




# mevo objects -----


#' Make a `mevo` object to store information needed for molecular evolution simulation.
#'
#'
#'
#' @param reference A \code{ref_genome} object from which you will generate variants.
#' @param sub Output from one of the \code{\link{sub_models}} functions that organizes
#'     information for the substitution models.
#'     See `?sub_models` for more information on these models and
#'     their required parameters.
#' @param ins Output from the \code{\link{indels}} function that specifies rates
#'     of insertions by length.
#'     Passing `NULL` to this argument results in no insertions.
#'     Defaults to `NULL`.
#' @param del Output from the \code{\link{indels}} function that specifies rates
#'     of deletions by length.
#'     Passing `NULL` to this argument results in no deletions.
#'     Defaults to `NULL`.
#' @param gamma_mats Output from the \code{\link{site_var}} function that specifies
#'     variability in mutation rates among sites (for both substitutions and indels).
#'     Passing `NULL` to this argument results in no variability among sites.
#'     Defaults to `NULL`.
#'
#' @return An object of class \code{\link{mevo}}.
#'
#' @noRd
#'
create_mevo <- function(reference,
                        sub,
                        ins,
                        del,
                        gamma_mats,
                        region_size) {

    if (!inherits(reference, "ref_genome")) {
        err_msg("create_variants", "reference", "a \"ref_genome\" object")
    }
    if (!is.null(sub) && !is_type(sub, "sub_model_info")) {
        err_msg("create_variants", "sub", "NULL or a \"sub_model_info\" object")
    }
    if (!is.null(ins) && !is_type(ins, "indel_rates")) {
        err_msg("create_variants", "ins", "NULL or a \"indel_rates\" object")
    }
    if (!is.null(del) && !is_type(del, "indel_rates")) {
        err_msg("create_variants", "del", "NULL or a \"indel_rates\" object")
    }
    if (!is.null(gamma_mats) && !is_type(gamma_mats, "site_var_mats")) {
        err_msg("create_variants", "gamma_mats", "NULL or a \"site_var_mats\" object")
    }
    if (!single_integer(region_size, 1)) {
        err_msg("create_variants", "region_size", "a single integer >= 1")
    }

    # `sub` must be provided if others are:
    if (is.null(sub) && (!is.null(ins) ||
                         !is.null(del) ||
                         !is.null(gamma_mats))) {
        stop("\nIn `create_variants`, if you want insertions, deletions, ",
             "or among-site variability, you must also provide substituion information ",
             "via one of the `sub_models` functions.", call. = FALSE)
    }

    # If no molecular evolution is provided, return NULL (only happens for VCF method)
    if (is.null(sub)) return(NULL)

    # Below will turn `NULL` into `numeric(0)` and
    # indel_rates object into simple numeric:
    ins <- as.numeric(ins)
    del <- as.numeric(del)


    # This results in no variability among sites and 1 Gamma region per sequence
    # (they'll get split later if desired):
    if (is.null(gamma_mats)) {
        seq_sizes <- reference$sizes()
        gamma_mats <- make_gamma_mats(seq_sizes, region_size_ = max(seq_sizes),
                                      shape = 0, invariant = 0)
        dim(gamma_mats) <- NULL # so it's just a list now
    }

    # -------+
    # Make final output object
    # -------+
    out <- mevo$new(sub,
                    ins,
                    del,
                    gamma_mats,
                    region_size)

    return(out)
}





#' Convert to a XPtr<MutationSampler> object
#'
#' @noRd
#'
mevo_obj_to_ptr <- function(mevo_obj) {

    sampler_ptr <- make_mutation_sampler_base(mevo_obj$Q,
                                              mevo_obj$pi_tcag,
                                              mevo_obj$insertion_rates,
                                              mevo_obj$deletion_rates,
                                              mevo_obj$region_size)

    return(sampler_ptr)
}



