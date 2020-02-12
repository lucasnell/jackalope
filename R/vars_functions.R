

#' Organize higher-level information for creating haplotypes.
#'
#'
#' The following functions organize information that gets passed to `create_haplotypes`
#' to generate haplotypes from a reference genome.
#' Each function represents a method of generation and starts with `"haps_"`.
#' The first three are phylogenomic methods, and all functions but `haps_vcf`
#' will use molecular evolution information when passed to `create_haplotypes`.
#'
#' \describe{
#'     \item{\code{\link{haps_theta}}}{Uses an estimate for theta, the population-scaled
#'         mutation rate, and a desired number of haplotypes.}
#'     \item{\code{\link{haps_phylo}}}{Uses phylogenetic tree(s) from `phylo`
#'         object(s) or NEWICK file(s), one tree per chromosome or one for all
#'         chromosomes.}
#'     \item{\code{\link{haps_gtrees}}}{Uses gene trees, either in the form of
#'         an object from the `scrm` or `coala` package or
#'         a file containing output in the style of the `ms` program.}
#'     \item{\code{\link{haps_ssites}}}{Uses matrices of segregating sites,
#'         either in the form of
#'         `scrm` or `coala` coalescent-simulator object(s), or
#'         a `ms`-style output file.}
#'     \item{\code{\link{haps_vcf}}}{Uses a haplotype call format (VCF) file that
#'         directly specifies haplotypes.}
#' }
#'
#'
#' @seealso \code{\link{create_haplotypes}}
#'
#' @name haps_functions
#'
NULL







# -------------*
#  Non-phylogenomic -----
# -------------*

#   __seg. sites -----

#' Process one segregating-sites matrix from a coalescent simulator with ms-style output.
#'
#' Used in `haps_ssites` below.
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



#' Organize information to create haplotypes using segregating sites matrices
#'
#'
#' This function organizes higher-level information for creating haplotypes from
#' matrices of segregating sites output from coalescent simulations.
#'
#'
#' For what the `seg_sites` field should look like in a list, see output from the
#' `scrm` or `coala` package.
#' (These packages are not required to be installed when installing
#' `jackalope`.)
#' If using either of these packages, I encourage you to cite them. For citation
#' information, see output from `citation("scrm")` or `citation("coala")`.
#'
#'
#' @param obj Object containing segregating sites information.
#'     This can be one of the following:
#'     (1) A single `list` with a `seg_sites` field inside. This field must
#'     contain a matrix for segregating sites for each chromosome.
#'     The matrix itself should contain the haplotype information, coded
#'     using 0s and 1s: 0s indicate the ancestral state and 1s indicate
#'     mutant.
#'     The matrix column names should indicate the positions of the polymorphisms on the
#'     chromosome.
#'     If positions are in the range `(0,1)`, they're assumed to come from an infinite-
#'     sites model and are relative positions.
#'     If positions are integers in the range `[0, chromosome length - 1]`
#'     or `[1, chromosome length]`, they're assumed to come from an finite-sites
#'     model and are absolute positions.
#'     Defaults to `NULL`.
#' @param fn A single string specifying the name of the file containing
#'     the `ms`-style coalescent output with segregating site info.
#'     Defaults to `NULL`.
#'
#'
#' @return A `haps_ssites_info` object containing information used in `create_haplotypes`
#'     to create variant haplotypes.
#'     This class is just a wrapper around a list of matrices of segregating site info,
#'     which you can view (but not change) using the object's `mats()` method.
#'
#' @export
#'
haps_ssites <- function(obj = NULL,
                        fn = NULL) {

    if (is.null(obj) && is.null(fn)) {
        stop("\nIn function `haps_ssites`, either argument `obj` or `fn` ",
             "must be provided.", call. = FALSE)
    }
    if (!is.null(obj) && !is.null(fn)) {
        stop("\nIn function `haps_ssites`, only one argument (`obj` or `fn`) ",
             "should be provided.", call. = FALSE)
    }

    sites_mats <- NULL
    if (!is.null(obj)) {

        # Check for coal_obj being a list and having a `seg_sites` field
        if (!inherits(obj, "list") || is.null(obj$seg_sites)) {
            err_msg("haps_ssites", "obj", "NULL or a list with a `seg_sites`",
                    "field present")
        }

        sites_mats <- lapply(obj$seg_sites, process_coal_obj_sites)

    } else {

        if (!is_type(fn, "character", 1)) {
            err_msg("haps_ssites", "fn", "NULL or a single string")
        }

        sites_mats <- coal_file_sites(fn)
        # Revert back to list (from arma::field which adds dims):
        dim(sites_mats) <- NULL
        n_cols <- sapply(sites_mats, ncol)
        if (any(n_cols < 2)) {
            stop("\nOne or more seg. sites matrices from a ms-style file output ",
                 "have no haplotype information specified.",
                 call. = FALSE)
        }

    }

    if (length(unique(sapply(sites_mats, ncol))) != 1) {
        stop("\nIn function `haps_ssites`, one or more of the segregating sites ",
             "matrices has a number of rows that differs from the rest.")
    }


    out <- haps_ssites_info$new(mats = sites_mats)

    return(out)

}






#   __vcf -----

#' Organize information to create haplotypes using a VCF file
#'
#' This function organizes higher-level information for creating haplotypes from
#' Variant Call Format (VCF) files.
#'
#'
#' @param fn A single string specifying the name of the VCF file
#' @param print_names Logical for whether to print all unique chromosome names from
#'     the VCF file when VCF chromosome names don't match those from the reference genome.
#'     This printing doesn't happen until this object is passed to `create_haplotypes`.
#'     This can be useful for troubleshooting.
#'     Defaults to `FALSE`.
#'
#' @export
#'
#' @return A `haps_vcf_info` object containing information used in `create_haplotypes`
#'     to create variant haplotypes.
#'     This class is just a wrapper around a list containing the arguments to this
#'     function, which you can view (but not change) using the object's `fn()` and
#'     `print_names()` methods.
#'
haps_vcf <- function(fn, print_names = FALSE) {

    if (!is_type(fn, "character", 1)) {
        err_msg("haps_vcf", "fn", "a single string")
    }
    if (!is_type(print_names, "logical", 1)) {
        err_msg("haps_vcf", "print_names", "a single logical")
    }

    fn <- path.expand(fn)

    if (!file.exists(fn)) {
        stop("\nFile ", fn, " doesn't exist.", call. = FALSE)
    }

    out <- haps_vcf_info$new(fn = fn, print_names = print_names)

    return(out)

}




# -------------*
#  Phylogenomic -----
# -------------*


# __phylo -----


#' Organize information to create haplotypes using phylogenetic tree(s)
#'
#' This function organizes higher-level information for creating haplotypes from
#' phylogenetic tree(s) output as `phylo` or `multiPhylo` objects
#' (both from the `ape` package) or NEWICK files.
#' Note that all phylogenetic trees must be rooted and binary.
#' If using this function, I encourage you to cite `ape`. For citation
#' information, see output from `citation("ape")`.
#'
#' See `?ape::write.tree` for writing phylogenies to an output file.
#'
#'
#' @param obj Object containing phylogenetic tree(s).
#'     This can be (1) a single `phylo` object
#'     that represents all chromosomes in the genome or
#'     (2) a `list` or `multiPhylo` object containing a `phylo` object for
#'     each reference chromosome.
#'     In the latter case, phylogenies will be assigned to chromosomes in the
#'     order provided.
#'     Defaults to `NULL`.
#' @param fn One or more string(s), each of which specifies the file name
#'     of a NEWICK file containing a phylogeny.
#'     If one name is provided, that phylogeny will be used for all chromosomes.
#'     If more than one is provided, there must be a phylogeny for each reference
#'     genome chromosome, and phylogenies will be assigned to chromosomes
#'     in the order provided.
#'     Defaults to `NULL`.
#'
#'
#' @return A `haps_phylo_info` object containing information used in `create_haplotypes`
#'     to create variant haplotypes.
#'     This class is just a wrapper around a list containing phylogenetic tree
#'     information for each reference chromosome, which you can view (but not change)
#'     using the object's `phylo()` method.
#'
#'
#' @export
#'
haps_phylo <- function(obj = NULL,
                       fn = NULL) {

    if (is.null(obj) && is.null(fn)) {
        stop("\nIn function `haps_phylo`, either argument `obj` or `fn` ",
             "must be provided.", call. = FALSE)
    }
    if (!is.null(obj) && !is.null(fn)) {
        stop("\nIn function `haps_phylo`, only one argument (`obj` or `fn`) ",
             "should be provided.", call. = FALSE)
    }

    phy <- NULL
    if (!is.null(obj)) {

        if ((!inherits(obj, "phylo") && !inherits(obj, "multiPhylo") &&
             !inherits(obj, "list")) ||
            (inherits(obj, "list") && !all(sapply(obj, inherits,
                                                  what = "phylo")))) {
            err_msg("haps_phylo", "obj",
                    "NULL or of class \"phylo\", \"multiPhylo\", or a list of",
                    "\"phylo\" objects")
        }

        phy <- obj
        # So all output is a list
        if (inherits(phy, "phylo")) phy <- list(phy)
        if (inherits(phy, "multiPhylo")) class(phy) <- "list"

    }
    if (!is.null(fn)) {
        if (!is_type(fn, "character")) {
            err_msg("haps_phylo", "fn", "NULL or a character vector")
        }
        phy <- lapply(fn, ape::read.tree)
    }

    out <- haps_phylo_info$new(phylo = phy)

    return(out)

}





# __theta -----

#' Organize information to create haplotypes using theta parameter
#'
#' This function organizes higher-level information for creating haplotypes from
#' the population-scaled mutation rate and a desired number of haplotypes.
#'
#'
#' @param theta Population-scaled mutation rate.
#' @param n_haps Number of desired haplotypes.
#'
#'
#' @return A `haps_theta_info` object containing information used in `create_haplotypes`
#'     to create variant haplotypes.
#'     This class is just a wrapper around a list containing the phylogenetic tree
#'     and `theta` parameter, which you can view (but not change) using the object's
#'     `phylo()` and `theta()` methods, respectively.
#'
#' @export
#'
haps_theta <- function(theta, n_haps) {

    if (!single_number(theta, 0)) {
        err_msg("haps_theta", "theta", "a single number >= 0")
    }
    if (!single_integer(n_haps, 2)) {
        err_msg("haps_theta", "n_haps", "a single integer >= 2")
    }


    # Generate random coalescent tree:
    phy <- ape::rcoal(n_haps)

    out <- haps_theta_info$new(phylo = phy, theta = theta)

    return(out)

}


# __gtrees -----

#' Organize information to create haplotypes using gene trees
#'
#' This function organizes higher-level information for creating haplotypes from
#' gene trees output from coalescent simulations.
#' Note that all gene trees must be rooted and binary.
#'
#' Using the `obj` argument is designed after the `trees` fields in the output from
#' the `scrm` and `coala` packages.
#' (These packages are not required to be installed when installing `jackalope`.)
#' To get gene trees, make sure to add `+ sumstat_trees()`
#' to the `coalmodel` for `coala`, or
#' make sure that `"-T"` is present in `args` for `scrm`.
#' If using either of these packages, I encourage you to cite them. For citation
#' information, see output from `citation("scrm")` or `citation("coala")`.
#'
#' If using an output file from a command-line program like `ms`/`msms`,
#' add the `-T` option.
#'
#'
#' @param obj Object containing gene trees.
#'     This can be one of the following:
#'     (1) A single `list` with a `trees` field inside. This field must
#'     contain a set of gene trees for each chromosome.
#'     (2) A list of lists, each sub-list containing a `trees` field of
#'     length 1. The top-level list must be of the same length as the
#'     number of chromosomes.
#'     Defaults to `NULL`.
#' @param fn A single string specifying the name of the file containing
#'     the `ms`-style coalescent output with gene trees.
#'     Defaults to `NULL`.
#'
#'
#' @return A `haps_gtrees_info` object containing information used in `create_haplotypes`
#'     to create variant haplotypes.
#'     This class is just a wrapper around a list of NEWICK tree strings, one for
#'     each gene tree, which you can view (but not change) using the object's
#'     `trees()` method.
#'
#' @export
#'
haps_gtrees <- function(obj = NULL,
                        fn = NULL) {

    if (is.null(obj) && is.null(fn)) {
        stop("\nIn function `haps_gtrees`, either argument `obj` or `fn` ",
             "must be provided.", call. = FALSE)
    }
    if (!is.null(obj) && !is.null(fn)) {
        stop("\nIn function `haps_gtrees`, only one argument (`obj` or `fn`) ",
             "should be provided.", call. = FALSE)
    }

    trees <- NULL

    if (!is.null(obj)) {

        # Check for `obj` being a list and either having a `trees` field or all its
        # items within having `trees` fields
        err <- FALSE
        nested <- FALSE
        if (!inherits(obj, "list")) {
            err <- TRUE
        } else if (is.null(obj$trees)) {
            if (any(sapply(obj, function(x) !inherits(x, "list") || is.null(x$trees)))) {
                err <- TRUE
            } else {
                nested <- TRUE
                if (any(sapply(obj, function(x) length(x$trees)) != 1)) err <- TRUE
                if (any(!sapply(obj, function(x) inherits(x$trees, "list")))) err <- TRUE
            }
        }
        if (err) {
            err_msg("haps_gtrees", "obj",
                    "(1) a list with a `trees` field present or",
                    "(2) a list of lists, each sub-list containing a `trees`",
                    "field of length 1")
        }

        if (nested) {
            trees <- lapply(obj, function(x) x$trees[[1]])
        } else trees <- obj$trees

    } else {

        if (!is_type(fn, "character", 1)) {
            err_msg("haps_gtrees", "fn", "NULL or a single string")
        }

        trees <- read_ms_trees_(fn)

        if (any(sapply(trees, length) == 0)) {
            stop("\nIn ms-style output file, one or more chromosomes have no trees.",
                 call. = FALSE)
        }

    }

    out <- haps_gtrees_info$new(trees = trees)

    return(out)

}



#' @describeIn haps_gtrees Write gene trees to ms-style output file.
#'
#' @param gtrees A `haps_gtrees_info` object output from `haps_gtrees`.
#' @param out_prefix Prefix for the output file of gene trees.
#'     The extension will be `.trees`.
#'
#' @export
#'
write_gtrees <- function(gtrees, out_prefix) {

    if (!inherits(gtrees, "haps_gtrees_info")) {
        err_msg("write_gtrees", "gtrees", "a `haps_gtrees_info` object")
    }
    if (!is_type(out_prefix, "character", 1L)) {
        err_msg("write_gtrees", "out_prefix", "a single string")
    }

    out_fn <- paste0(path.expand(out_prefix), ".trees")

    # Trees separated by empty line, then one with just "//"
    # This is to emulate ms-style output
    out_str <- paste0("//\n", paste(lapply(gtrees$trees(), paste, collapse = "\n"),
                                    collapse = "\n\n//\n"))

    write(out_str, out_fn)

    invisible(NULL)

}


