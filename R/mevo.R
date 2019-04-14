
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

#   SEG. SITES METHOD -----

# ====================================================================================`
# ====================================================================================`


#' Process one segregating-sites matrix from a coalescent simulator with ms-style output.
#'
#' Used in `vars_ssites` below.
#'
#' @param mat The matrix to process.
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


#' Check validity of position columns in segregating-sites matrices.
#'
#' Used in `to_var_set.vars_ssites_info` below.
#'
#' @noRd
#'
fill_coal_mat_pos <- function(sites_mats, seq_sizes) {

    if (length(sites_mats) != length(seq_sizes)) {
        stop("\nIn function `vars_ssites`, there must be exactly one segregating sites ",
             "matrix for each reference genome sequence. ",
             "It appears you need to re-run `vars_ssites` before attempting to ",
             "run `create_variants` again.")
    }

    for (i in 1:length(sites_mats)) {
        if (nrow(sites_mats[[i]]) == 0) next;
        pos <- sites_mats[[i]][,1]
        if (all(pos < 1 && pos > 0)) {
            # Converting to integer positions (0-based):
            pos  <- as.integer(pos * seq_sizes[i]);
        } else {
            all_ints <- all(pos %% 1 == 0)
            if (all_ints && all(pos < seq_sizes[i] && pos >= 0)) {
                # Keeping them in 0-based indices:
                pos = as.integer(pos);
            } else if (all_ints && all(pos <= seq_sizes[i] && pos >= 1)) {
                # Converting to 0-based indices:
                pos = as.integer(pos) - 1;
            } else {
                stop("\nPositions in one or more segregating-sites matrices ",
                     "are not obviously from either a finite- or infinite-sites model. ",
                     "The former should have integer positions in the range ",
                     "[0, sequence length - 1] or [1, sequence length], ",
                     "the latter numeric in (0,1).",
                     "It appears you need to re-run `vars_ssites` before attempting to ",
                     "run `create_variants` again.")
            }
        }
        sites_mats[[i]][,1] <- pos

    }

    return(sites_mats)

}




#' Create necessary information to create variants using segregating sites matrices
#'
#'
#' This function organizes higher-level information for creating variants from
#' matrices of segregating sites output from coalescent simulations.
#'
#'
#' For what the `seg_sites` field should look like in a list, see output from the
#' `scrm` or `coala` package.
#' (These packages are not required to be installed when installing
#' `jackalope`.)
#'
#' @param obj Object containing segregating sites information.
#'     This can be one of the following:
#'     (1) A single `list` with a `seg_sites` field inside. This field must
#'     contain a matrix for segregating sites for each sequence.
#'     The matrix itself should contain the haplotype information, coded
#'     using 0s and 1s: 0s indicate the ancestral state and 1s indicate
#'     mutant.
#'     The matrix column names should indicate the positions of the polymorphisms on the
#'     chromosome.
#'     If positions are in the range `(0,1)`, they're assumed to come from an infinite-
#'     sites model and are relative positions.
#'     If positions are integers in the range `[0, sequence length - 1]`
#'     or `[1, sequence length]`, they're assumed to come from an finite-sites
#'     model and are absolute positions.
#'     Defaults to `NULL`.
#' @param fn A single string specifying the name of the file containing
#'     the `ms`-style coalescent output with segregating site info.
#'     Defaults to `NULL`.
#'
#'
#' @return A `vars_ssites_info` object containing information used in `create_variants`
#'     to create haploid variants.
#'     This class is just a wrapper around a list of matrices of segregating site info.
#'
#' @export
#'
vars_ssites <- function(obj = NULL,
                        fn = NULL) {

    if (is.null(obj) && is.null(fn)) {
        stop("\nIn function `vars_ssites`, either argument `obj` or `fn` ",
             "must be provided.", call. = FALSE)
    }
    if (!is.null(obj) && !is.null(fn)) {
        stop("\nIn function `vars_ssites`, only one argument (`obj` or `fn`) ",
             "should be provided.", call. = FALSE)
    }

    sites_mats <- NULL
    if (!is.null(obj)) {

        # Check for coal_obj being a list and having a `seg_sites` field
        if (!inherits(obj, "list") || is.null(obj$seg_sites)) {
            err_msg("vars_ssites", "obj", "a list with a `seg_sites` field present")
        }

        sites_mats <- lapply(obj$seg_sites, process_coal_obj_sites)

    } else {

        if (!is_type(fn, "character", 1)) err_msg("vars_ssites", "fn", "a single string")

        sites_mats <- coal_file_sites(fn)
        # Revert back to list (from arma::field which adds dims):
        dim(sites_mats) <- NULL
        n_cols <- sapply(sites_mats, ncol)
        if (any(n_cols < 2)) {
            stop("\nOne or more seg. sites matrices from a ms-style file output ",
                 "have no variant information specified.",
                 call. = FALSE)
        }

    }

    if (length(unique(sapply(sites_mats, ncol))) != 1) {
        stop("\nIn function `vars_ssites`, one or more of the segregating sites ",
             "matrices has a number of rows that differs from the rest.")
    }


    out <- list(mats = sites_mats)
    class(out) <- "vars_ssites_info"

    return(out)

}

#' Print output from `vars_ssites`
#'
#' @noRd
#' @export
#'
print.vars_ssites_info <- function(x, digits = max(3, getOption("digits") - 3), ...) {

    cat("< Seg. site variant-creation info >\n")
    cat(sprintf("  * Number of variants: %i\n", ncol(x$mats[[1]]) - 1))
    cat(sprintf("  * Number of sites: %s\n", format(as.integer(sum(sapply(x$mats, nrow))),
                                                    big.mark = ",")))
    invisible(NULL)

}


#' Create variants from segregating-site info from coalescent simulations.
#'
#'
#' @noRd
#'
to_var_set.vars_ssites_info <- function(x, reference, mevo_obj,
                                        n_threads, show_progress, ...) {


    seq_sizes <- reference$sizes()

    # Fill and check the position column in `x$mats`
    x$mats <- fill_coal_mat_pos(x$mats, seq_sizes)

    variants_ptr <- add_coal_sites_cpp(reference$genome,
                                       x$mats,
                                       mevo_obj$Q,
                                       mevo_obj$pi_tcag,
                                       mevo_obj$insertion_rates,
                                       mevo_obj$deletion_rates,
                                       n_threads,
                                       show_progress)

    return(variants_ptr)

}







# ====================================================================================`
# ====================================================================================`

#   VCF METHOD -----

# ====================================================================================`
# ====================================================================================`

#' Create necessary information to create variants using a VCF file
#'
#' This function organizes higher-level information for creating variants from
#' Variant Call Format (VCF) files.
#'
#' This function won't work if the package `vcfR` isn't installed.
#'
#' @param fn A single string specifying the name of the VCF file
#' @param print_names Logical for whether to print all unique sequence names from
#'     the VCF file when VCF sequence names don't match those from the reference genome.
#'     This printing doesn't happen until this object is passed to `create_variants`.
#'     This can be useful for troubleshooting.
#'     Defaults to `FALSE`.
#' @param ... Arguments to pass to `vcfR::read.vcfR`, excluding the `file` argument
#'     that will be overridden with the `fn` argument to this function.
#'
#' @export
#'
#' @return A `vars_vcf_info` object containing information used in `create_variants`
#'     to create haploid variants.
#'     This class is just a wrapper around a list containing relevant output from
#'     `vcfR::read.vcfR`:
#'     haplotypes, reference sequences, positions, sequence names, and variant names.
#'
vars_vcf <- function(fn, print_names = FALSE, ...) {

    if (!requireNamespace("vcfR", quietly = TRUE)) {
        stop("\nPackage \"vcfR\" is needed for reading VCF files. ",
             "Please install it.",
             call. = FALSE)
    }

    if (!is_type(fn, "character", 1)) {
        err_msg("vars_vcf", "fn", "a single string")
    }

    other_args <- list(...)
    if (is.null(other_args$verbose)) other_args$verbose <- FALSE
    method_info$file <- fn

    if (!all(names(other_args) %in% names(formals(vcfR::read.vcfR)))) {
        bad_names <- names(other_args)[!names(other_args) %in%
                                           names(formals(vcfR::read.vcfR))]
        stop("\nIn function `vars_vcf` in jackalope, the following extra ",
             "arguments provided don't match any arguments in `vcfR::read.vcfR`: ",
             paste(bad_names, collapse = ", "), ".", call. = FALSE)
    }

    vcf <- do.call(vcfR::read.vcfR, read_args)

    chrom <- vcf@fix[,"CHROM"]
    pos <- as.integer(vcf@fix[,"POS"]) - 1  # -1 is to convert to C++ indices
    ref_seq <- vcf@fix[,"REF"]
    alts <- strsplit(vcf@fix[,"ALT"], ",")

    if (length(chrom) != length(pos) | length(chrom) != length(ref_seq)) {
        stop("\nVCF not parsing correctly. ",
             "Vectors of chromosomes, positions, and reference-sequences aren't ",
             "all the same length.",
             call. = FALSE)
    }

    haps <- vcfR::extract.haps(vcf, unphased_as_NA = FALSE, verbose = FALSE)

    if (length(chrom) != length(pos) | length(chrom) != length(ref_seq)) {
        stop("\nVCF not parsing correctly. ",
             "Vectors of chromosomes, positions, and reference-sequences aren't ",
             "all the same length.",
             call. = FALSE)
    }
    if (nrow(haps) != length(pos)) {
        stop("\nVCF not parsing correctly. ",
             "Number of haplotypes doesn't match with number of positions.",
             call. = FALSE)
    }


    # I'm assuming NAs mean no mutation
    haps[is.na(haps)] <- ""

    var_names <- colnames(haps)
    colnames(haps) <- NULL
    rownames(haps) <- NULL

    # Split into list for easier processing in `read_vcfr`
    haps <- split(haps, row(haps))

    # We treat things differently if vcfR has output numbers rather than nucleotides.
    # (The below line should be TRUE when it outputs numbers.)
    if (any(!is.na(suppressWarnings(as.integer(do.call(c, haps)))))) {
        # Change string integers to actual genotypes:
        haps <-
            lapply(1:length(haps),
                   function(i) {
                       as.character(sapply(haps[[i]],
                                           function(j) {
                                               ifelse(j == "" | j == "0", "",
                                                      alts[[i]][as.integer(j)])
                                           }))
                   })
    }

    out <- list(haps = haps, pos = pos, chrom = chrom, var_names = var_names,
                ref_seq = ref_seq)

    class(out) <- "vars_vcf_info"

    return(out)

}


#' Print output from `vars_vcf`
#'
#'
#' @noRd
#' @export
#'
print.vars_vcf_info <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("< VCF variant-creation info >\n")
    cat(sprintf("  * Number of variants: %i\n", ncol(x$haps[[1]])))
    cat(sprintf("  * Number of sequences: %i\n", length(unique(x$chrom))))
    cat(sprintf("  * Number of sites: %s\n", format(as.integer(length(x$pos)),
                                                    big.mark = ",")))
    invisible(NULL)
}


#' Create variants from VCF file
#'
#'
#' @noRd
#'
to_var_set.vars_vcf_info <- function(x, reference, mevo_obj, n_threads, show_progress) {

    seq_names <- view_ref_genome_seq_names(reference$genome)
    unq_chrom <- unique(x$chrom)

    if (!all(unq_chrom %in% seq_names)) {
        if (print_names) print(unq_chrom)
        stop("\nSequence name(s) in VCF file don't match those in the ",
             "`ref_genome` object. ",
             "It's probably easiest to manually change the `ref_genome` object ",
             "(using `$set_names()` method) to have the same names as the VCF file. ",
             "Re-run this function with `print_names = TRUE` to see the VCF-file names.",
             call. = FALSE)
    }

    # Converts items in `chrom` to 0-based indices of sequences in ref. genome
    chrom_inds <- match(x$chrom, seq_names) - 1


    variants_ptr <- read_vcfr(reference$genome, x$var_names,
                              x$haps, chrom_inds, x$pos, x$ref_seq)

    return(variants_ptr)

}
