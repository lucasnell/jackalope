
# ====================================================================================`
# ====================================================================================`

#   MEVO OBJECT -----

# ====================================================================================`
# ====================================================================================`

#' Print method for sub_model_info objects.
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
    if (!single_integer(max_length, 1)) {
        err_msg("indels", "max_length", "a single integer >= 1")
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
#'
#' @param shape Shape parameter for the Gamma distribution that generates gamma distances,
#'     The variance of the distribution is `1 / shape`, and its mean is fixed to 1.
#'     Defaults to `NULL`.
#' @param region_size Size of regions to break the genome into,
#'     where all sites within a region have the same gamma distance.
#'     Defaults to `NULL`.
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
#'                            rgamma(10, 0.9))))
#'
site_var <- function(reference,
                     shape = NULL,
                     region_size = NULL,
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
        mats <- make_gamma_mats(seq_sizes, gamma_size_ = region_size, shape = shape)
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


#' Print method for site_var_mats objects.
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
#' @noRd
#'
create_mevo <- function(reference,
                        sub,
                        ins,
                        del,
                        gamma_mats,
                        chunk_size) {

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
    if (!single_integer(chunk_size, 0)) {
        err_msg("create_variants", "chunk_size", "an integer >= 0")
    }

    # `sub` must be provided if others are:
    if (is.null(sub) && (!is.null(ins) ||
                         !is.null(del) ||
                         !is.null(gamma_mats))) {
        stop("\nIn `create_variants`, if you want insertions, deletions, ",
             "or among-site variability, you must also provide substituion information ",
             "via one of the `sub_models` functions.", call. = FALSE)
    }

    # If no molecular evolution is provided, return NULL
    if (is.null(sub)) return(NULL)

    # Below will turn `NULL` into `numeric(0)` and
    # indel_rates object into simple numeric:
    ins <- as.numeric(ins)
    del <- as.numeric(del)


    # -------+
    # Process info for mutation-rate variability among sites and write to BED
    # file if desired
    # -------+
    if (is.null(gamma_mats)) {
        # This results in no variability among sites:
        gamma_mats <- make_gamma_mats(seq_sizes, gamma_size_ = 0, shape = 1)
        dim(gamma_mats) <- NULL # so it's just a list now
    }

    # -------+
    # Make final output object
    # -------+
    out <- mevo$new(sub,
                    ins,
                    del,
                    gamma_mats,
                    chunk_size)

    return(out)
}




# ====================================================================================`
# ====================================================================================`

#   COAL SITES -----

# ====================================================================================`
# ====================================================================================`



#' Process one segregating-sites matrix from a coalescent simulator with ms-style output.
#'
#' @param mat The matrix to process.
#' @param seq_size The number of bp in the sequence associated with the input string.
#'
#' @noRd
#'
process_coal_obj_sites <- function(mat) {

    err_msg_ <- paste("\nPositions in one or more segregating-sites matrices %s.",
                      "They are derived from column names, so check those for nonsense.")

    # Dealing with coala objects:
    if (inherits(mat, "segsites")) mat <- mat$snps

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
                "`names` field")
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
            stop("\nYour `method_info` argument to `create_variants` is causing ",
                 "problems: ",
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
read_coal_sites <- function(method_info, reference, mevo_obj, n_threads, show_progress) {

    var_names <- character(0)
    if (inherits(method_info, "list") && !is.null(method_info$names) &&
        !is.null(method_info$info)) {
        var_names <- method_info$names
        method_info <- method_info$info
    }

    if (inherits(method_info, "list")) {
        sites_mats <- coal_obj_sites(method_info)
    } else if (is_type(method_info, "character", 1)) {
        sites_mats <- coal_file_sites(method_info)
        # Revert back to list (from arma::field which adds dims)
        dim(sites_mats) <- NULL
        n_cols <- sapply(sites_mats, ncol)
        if (any(n_cols < 2)) {
            stop("\nOne or more seg. sites matrices from a ms-style file output ",
                 "have no variant information specified.",
                 call. = FALSE)
        }
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
                                       n_threads,
                                       show_progress)

    return(variants_ptr)

}


