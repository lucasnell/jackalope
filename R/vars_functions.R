

#' Organize higher-level information for creating variants.
#'
#'
#' The following functions organize information that gets passed to `create_variants`
#' to generate variants from a reference genome.
#' Each function represents a method of generation and starts with `"vars_"`.
#' The first three are phylogenomic methods, and all functions but `vars_vcf`
#' will use molecular evolution information when passed to `create_variants`.
#'
#' \describe{
#'     \item{\code{\link{vars_theta}}}{Uses an estimate for theta, the population-scaled
#'         mutation rate, and a desired number of variants.}
#'     \item{\code{\link{vars_phylo}}}{Uses phylogenetic tree(s) from `phylo`
#'         object(s) or NEWICK file(s), one tree per chromosome or one for all
#'         chromosomes.}
#'     \item{\code{\link{vars_gtrees}}}{Uses gene trees, either in the form of
#'         an object from the `scrm` or `coala` package or
#'         a file containing output in the style of the `ms` program.}
#'     \item{\code{\link{vars_ssites}}}{Uses matrices of segregating sites,
#'         either in the form of
#'         `scrm` or `coala` coalescent-simulator object(s), or
#'         (2) a `ms`-style output file.}
#'     \item{\code{\link{vars_vcf}}}{Uses a variant call format (VCF) file that
#'         directly specifies variants.
#'         This method does not work if the `vcfR` package isn't installed.}
#' }
#'
#'
#' @seealso \code{\link{create_variants}}
#'
#' @name vars_functions
#'
NULL





# ====================================================================================`
# ====================================================================================`

# * MAIN * -----

# ====================================================================================`
# ====================================================================================`


# -------------*
#  Non-phylogenomic -----
# -------------*

#   __seg. sites -----

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
            err_msg("vars_ssites", "obj", "NULL or a list with a `seg_sites`",
                    "field present")
        }

        sites_mats <- lapply(obj$seg_sites, process_coal_obj_sites)

    } else {

        if (!is_type(fn, "character", 1)) {
            err_msg("vars_ssites", "fn", "NULL or a single string")
        }

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






#   __vcf -----

#' Create necessary information to create variants using a VCF file
#'
#' This function organizes higher-level information for creating variants from
#' Variant Call Format (VCF) files.
#'
#' This function won't work if the package `vcfR` isn't installed.
#' If using this function, I encourage you to cite `vcfR`. For citation
#' information, see output from `citation("vcfR")`.
#'
#' @param fn A single string specifying the name of the VCF file
#' @param print_names Logical for whether to print all unique chromosome names from
#'     the VCF file when VCF chromosome names don't match those from the reference genome.
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
#'     haplotypes, reference chromosomes, positions, chromosome names, and variant names.
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

    read_args <- list(...)
    if (is.null(read_args$verbose)) read_args$verbose <- FALSE
    read_args$file <- fn

    if (!all(names(read_args) %in% names(formals(vcfR::read.vcfR)))) {
        bad_names <- names(read_args)[!names(read_args) %in%
                                           names(formals(vcfR::read.vcfR))]
        stop("\nIn function `vars_vcf` in jackalope, the following extra ",
             "arguments provided don't match any arguments in `vcfR::read.vcfR`: ",
             paste(bad_names, collapse = ", "), ".", call. = FALSE)
    }

    vcf <- do.call(vcfR::read.vcfR, read_args)

    chrom <- vcf@fix[,"CHROM"]
    pos <- as.integer(vcf@fix[,"POS"]) - 1  # -1 is to convert to C++ indices
    ref_chrom <- vcf@fix[,"REF"]
    alts <- strsplit(vcf@fix[,"ALT"], ",")

    if (length(chrom) != length(pos) | length(chrom) != length(ref_chrom)) {
        stop("\nVCF not parsing correctly. ",
             "Vectors of chromosomes, positions, and reference-chromosomes aren't ",
             "all the same length.",
             call. = FALSE)
    }

    haps <- vcfR::extract.haps(vcf, unphased_as_NA = FALSE, verbose = FALSE)

    if (length(chrom) != length(pos) | length(chrom) != length(ref_chrom)) {
        stop("\nVCF not parsing correctly. ",
             "Vectors of chromosomes, positions, and reference-chromosomes aren't ",
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
                ref_chrom = ref_chrom, print_names = print_names)

    class(out) <- "vars_vcf_info"

    return(out)

}




# -------------*
#  Phylogenomic -----
# -------------*


# __phylo -----


#' Create necessary information to create variants using phylogenetic tree(s)
#'
#' This function organizes higher-level information for creating variants from
#' phylogenetic tree(s) output as `phylo` or `multiPhylo` objects
#' (both from the `ape` package) or NEWICK files.
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
#' @return A `vars_phylo_info` object containing information used in `create_variants`
#'     to create haploid variants.
#'     This class is just a wrapper around a list containing phylogenetic tree
#'     information for each reference chromosome.
#'
#'
#' @export
#'
vars_phylo <- function(obj = NULL,
                       fn = NULL) {

    if (is.null(obj) && is.null(fn)) {
        stop("\nIn function `vars_phylo`, either argument `obj` or `fn` ",
             "must be provided.", call. = FALSE)
    }
    if (!is.null(obj) && !is.null(fn)) {
        stop("\nIn function `vars_phylo`, only one argument (`obj` or `fn`) ",
             "should be provided.", call. = FALSE)
    }

    phy <- NULL
    if (!is.null(obj)) {

        if ((!inherits(obj, "phylo") && !inherits(obj, "multiPhylo") &&
             !inherits(obj, "list")) ||
            (inherits(obj, "list") && !all(sapply(obj, inherits,
                                                  what = "phylo")))) {
            err_msg("vars_phylo", "obj",
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
            err_msg("vars_phylo", "fn", "NULL or a character vector")
        }
        phy <- lapply(fn, ape::read.tree)
    }

    out <- list(phylo = phy)
    class(out) <- "vars_phylo_info"

    return(out)

}





# __theta -----

#' Create necessary information to create variants using theta parameter
#'
#' This function organizes higher-level information for creating variants from
#' the population-scaled mutation rate and a desired number of variants.
#'
#'
#' @param theta Population-scaled mutation rate.
#' @param n_vars Mumber of desired variants.
#'
#'
#' @return A `vars_theta_info` object containing information used in `create_variants`
#'     to create haploid variants.
#'     This class is just a wrapper around a list containing the phylogenetic tree
#'     and `theta` parameter.
#'
#' @export
#'
vars_theta <- function(theta, n_vars) {

    if (!single_number(theta, 0)) {
        err_msg("vars_theta", "theta", "a single number >= 0")
    }
    if (!single_integer(n_vars, 2)) {
        err_msg("vars_theta", "n_vars", "a single integer >= 2")
    }


    # Generate random coalescent tree:
    phy <- ape::rcoal(n_vars)

    out <- list(phylo = phy, theta = theta)
    class(out) <- "vars_theta_info"

    return(out)

}


# __gtrees -----

#' Create necessary information to create variants using gene trees
#'
#' This function organizes higher-level information for creating variants from
#' gene trees output from coalescent simulations.
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
#' @return A `vars_gtrees_info` object containing information used in `create_variants`
#'     to create haploid variants.
#'     This class is just a wrapper around a list of NEWICK tree strings, one for
#'     each gene tree.
#'
#' @export
#'
vars_gtrees <- function(obj = NULL,
                        fn = NULL) {

    if (is.null(obj) && is.null(fn)) {
        stop("\nIn function `vars_gtrees`, either argument `obj` or `fn` ",
             "must be provided.", call. = FALSE)
    }
    if (!is.null(obj) && !is.null(fn)) {
        stop("\nIn function `vars_gtrees`, only one argument (`obj` or `fn`) ",
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
            err_msg("vars_gtrees", "obj",
                    "(1) a list with a `trees` field present or",
                    "(2) a list of lists, each sub-list containing a `trees`",
                    "field of length 1")
        }

        if (nested) {
            trees <- lapply(obj, function(x) x$trees[[1]])
        } else trees <- obj$trees

    } else {

        if (!is_type(fn, "character", 1)) {
            err_msg("vars_gtrees", "fn", "NULL or a single string")
        }

        trees <- read_ms_trees_(fn)

        if (any(sapply(trees, length) == 0)) {
            stop("\nIn ms-style output file, one or more chromosomes have no trees.",
                 call. = FALSE)
        }

    }

    out <- list(trees = trees)
    class(out) <- "vars_gtrees_info"

    return(out)

}



#' @describeIn vars_gtrees Write gene trees to ms-style output file.
#'
#' @param gtrees A `vars_gtrees_info` object output from `vars_gtrees`.
#' @param out_prefix Prefix for the output file of gene trees.
#'     The extension will be `.trees`.
#'
#' @export
#'
write_gtrees <- function(gtrees, out_prefix) {

    if (!inherits(gtrees, "vars_gtrees_info")) {
        err_msg("write_gtrees", "gtrees", "a `vars_gtrees_info` object")
    }
    if (!is_type(out_prefix, "character", 1L)) {
        err_msg("write_gtrees", "out_prefix", "a single string")
    }

    out_fn <- paste0(path.expand(out_prefix), ".trees")

    # Trees separated by empty line, then one with just "//"
    # This is to emulate ms-style output
    out_str <- paste0("//\n", paste(lapply(gtrees$trees, paste, collapse = "\n"),
                                    collapse = "\n\n//\n"))

    write(out_str, out_fn)

    invisible(NULL)

}



# ====================================================================================`
# ====================================================================================`

# * HELPERS * -----

# ====================================================================================`
# ====================================================================================`


# -------------*
#  Non-phylogenomic -----
# -------------*

#   __seg. sites -----


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
fill_coal_mat_pos <- function(sites_mats, chrom_sizes) {

    if (length(sites_mats) != length(chrom_sizes)) {
        stop("\nIn function `vars_ssites`, there must be exactly one segregating sites ",
             "matrix for each reference genome chromosome. ",
             "It appears you need to re-run `vars_ssites` before attempting to ",
             "run `create_variants` again.")
    }

    for (i in 1:length(sites_mats)) {
        if (nrow(sites_mats[[i]]) == 0) next;
        pos <- sites_mats[[i]][,1]
        if (all(pos < 1 & pos > 0)) {
            # Converting to integer positions (0-based):
            pos  <- as.integer(pos * chrom_sizes[i]);
        } else {
            all_ints <- all(pos %% 1 == 0)
            if (all_ints && all(pos < chrom_sizes[i] & pos >= 0)) {
                # Keeping them in 0-based indices:
                pos = as.integer(pos);
            } else if (all_ints && all(pos <= chrom_sizes[i] & pos >= 1)) {
                # Converting to 0-based indices:
                pos = as.integer(pos) - 1;
            } else {
                stop("\nPositions in one or more segregating-sites matrices ",
                     "are not obviously from either a finite- or infinite-sites model. ",
                     "The former should have integer positions in the range ",
                     "[0, chromosome length - 1] or [1, chromosome length], ",
                     "the latter numeric in (0,1).",
                     "It appears you need to re-run `vars_ssites` before attempting to ",
                     "run `create_variants` again.")
            }
        }
        sites_mats[[i]][,1] <- pos

    }

    return(sites_mats)

}




# -------------*
#  Phylogenomic -----
# -------------*


# These first two are helpers for multiple phylogenomic methods


#' Go from pointer to trees info to a pointer to a VarSet object
#'
#' Used below in `to_var_set` methods
#'
#' @noRd
#'
trees_to_var_set <- function(phylo_info_ptr, reference, sub, ins, del, epsilon,
                             n_threads, show_progress) {

    # Make sampler_base_ptr
    sampler_base_ptr <- make_mutation_sampler_base(sub$Q(),
                                                   sub$U(),
                                                   sub$Ui(),
                                                   sub$L(),
                                                   sub$invariant(),
                                                   ins$rates(),
                                                   del$rates(),
                                                   epsilon,
                                                   sub$pi_tcag())

    # Make variants pointer:
    variants_ptr <- evolve_chroms(reference$ptr(),
                                  sampler_base_ptr,
                                  phylo_info_ptr,
                                  n_threads,
                                  show_progress)

    return(variants_ptr)

}




#' Read info from a `phylo` object.
#'
#' @return An external pointer to the phylogenetic info needed to do the simulations.
#'
#' @noRd
#'
phylo_to_ptr <- function(phy, n_chroms) {

    if ((!inherits(phy, "phylo") && !inherits(phy, "multiPhylo") &&
         !inherits(phy, "list")) ||
        (inherits(phy, "list") && !all(sapply(phy, inherits, what = "phylo")))) {
        stop("\nThe `phy` argument to the internal function `phylo_to_ptr` should ",
             "only ever be an object of class \"phylo\", \"multiPhylo\", or a ",
             "list of \"phylo\" objects.")
    }

    if ((inherits(phy, "multiPhylo") || inherits(phy, "list")) &&
        length(phy) != n_chroms) {
        stop("\nThe `phy` argument to the internal function `phylo_to_ptr` should ",
             "have a length of 1 or equal to the number of chromosomes if it's of class",
             "\"multiPhylo\" or list.")
    }
    # So all inputs are lists of the proper length:
    if (inherits(phy, "phylo")) phy <- rep(list(phy), n_chroms)
    if (inherits(phy, "multiPhylo")) class(phy) <- "list"

    phylo_info <- lapply(phy,
                         function(p) {
                             p <- ape::reorder.phylo(p, order = "cladewise")
                             labels <- paste(p$tip.label)# <-- making sure they're strings
                             branch_lens <- p$edge.length
                             edges <- p$edge
                             phy_info <- list(branch_lens = branch_lens, edges = edges,
                                              labels = labels, start = 0, end = 0)
                             return(list(phy_info))
                         })

    trees_ptr <- phylo_info_to_trees(phylo_info)

    return(trees_ptr)
}



# __gtrees -----


#' Process one gene-tree string from a coalescent simulator with ms-style output.
#'
#' Used in `gtrees_to_ptr`.
#'
#' @param str The string to process.
#' @param chrom_size The number of bp in the chromosome associated with the input string.
#'
#' @noRd
#'
process_coal_tree_string <- function(str, chrom_size) {

    if (length(str) > 1) {
        if (!all(grepl("^\\[", str))) {
            stop("\nA coalescent string appears to include ",
                 "recombination but does not include sizes for each region.",
                 call. = FALSE)
        }
        sizes_ <- as.numeric(sapply(str, function(x) strsplit(x, "\\[|\\]")[[1]][2]))
        # If they're <= 1, then they're not # bp, they're proportion of chromosome
        if (all(sizes_ <= 1) & chrom_size > 1) {
            sizes_ <- sizes_ / sum(sizes_)
            sizes_ <- round(sizes_ * chrom_size, 0)
            # Remove any zero sizes:
            sizes_ <- sizes_[sizes_ > 0]
            # If there's nothing left, just make it of length 1:
            if (length(sizes_) == 0) {
                sizes_ <- chrom_size
                # If it doesn't round quite right, then randomly add/subtract:
            } else if (sum(sizes_) != chrom_size) {
                inds <- sample.int(length(sizes_), abs(chrom_size - sum(sizes_)))
                sizes_[inds] <- sizes_[inds] + sign(chrom_size - sum(sizes_))
            }
        } else if (sum(sizes_) != chrom_size) {
            stop("\nA coalescent string appears to include ",
                 "recombination but the combined sizes of all regions don't match ",
                 "the size of the chromosome.", call. = FALSE)
        }
    } else {
        sizes_ <- 1
    }
    ends <- cumsum(sizes_) - 1
    starts <- c(0, utils::head(ends, -1) + 1)

    phylo_ <- ape::read.tree(text = str)
    # If no recombination (so only one phylo per chromosome),
    # then we need to make sure that it has the same nestedness as if there were
    # >1 phylo objects:
    if (inherits(phylo_, "phylo")) {
        phylo_ <- list(phylo_)
    }

    out <- rep(list(NA), length(phylo_))

    for (i in 1:length(phylo_)) {
        phy <- ape::reorder.phylo(phylo_[[i]], order = "cladewise")
        labels <- paste(phy$tip.label) # used paste to make sure they're characters
        branch_lens <- phy$edge.length
        edges <- phy$edge
        out[[i]] <- list(branch_lens = branch_lens, edges = edges, labels = labels,
                         start = starts[i], end = ends[i])
    }

    return(out)
}


#' Create pointer to trees C++ class from gene-tree info.
#'
#' Used in `to_var_set.vars_gtrees_info`.
#'
#' @return An XPtr to the info needed from the gene trees to do the simulations.
#'
#' @noRd
#'
#'
gtrees_to_ptr <- function(trees, reference) {

    chrom_sizes <- reference$sizes()

    if (length(trees) != length(chrom_sizes)) {
        stop("\nIn function `vars_gtrees`, there must be a set of gene trees ",
             "for each reference genome chromosome. ",
             "It appears you need to re-run `vars_gtrees` before attempting to ",
             "run `create_variants` again.")
    }

    phylo_info <- mapply(process_coal_tree_string, trees, chrom_sizes,
                         SIMPLIFY = FALSE, USE.NAMES = FALSE)

    unq_n_tips <- lapply(phylo_info,
                         function(x) sapply(x, function(xx) length(xx$labels)))
    unq_n_tips <- unique(do.call(c, unq_n_tips))
    if (length(unq_n_tips) > 1) {
        stop("\nIn function `vars_gtrees`, all gene trees must have the same ",
             "number of tips. ",
             "It appears you need to re-run `vars_gtrees` before attempting to ",
             "run `create_variants` again.")
    }
    unq_tips_names <- sapply(phylo_info,
                             function(x) {
                                 tips_ <- do.call(c, lapply(x, function(xx) xx$labels))
                                 paste(sort(unique(tips_)), collapse = "___")
                             })
    if (length(unique(unq_tips_names)) > 1) {
        stop("\nIn function `vars_gtrees`, all gene trees must have the same ",
             "tip names. ",
             "It appears you need to re-run `vars_gtrees` before attempting to ",
             "run `create_variants` again.")
    }

    # Making sure all labels are the same
    label_mat <- do.call(rbind,
                         lapply(phylo_info, function(x) {
                             t(sapply(x, function(xx) xx$labels))
                         }))
    label_mat <- t(apply(label_mat, 1, sort))
    for (i in 2:nrow(label_mat)) {
        if (any(label_mat[1,] != label_mat[i,])) {
            stop("\nIn function `vars_gtrees`, all gene trees must have the same ",
                 "tip names. ",
                 "It appears you need to re-run `vars_gtrees` before attempting to ",
                 "run `create_variants` again.")
        }
    }

    trees_ptr <- phylo_info_to_trees(phylo_info)

    return(trees_ptr)
}






# ====================================================================================`
# ====================================================================================`

# * PRINT * -----

# ====================================================================================`
# ====================================================================================`

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


#' Print output from `vars_vcf`
#'
#'
#' @noRd
#' @export
#'
print.vars_vcf_info <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("< VCF variant-creation info >\n")
    cat(sprintf("  * Number of variants: %i\n", ncol(x$haps[[1]])))
    cat(sprintf("  * Number of chromosomes: %i\n", length(unique(x$chrom))))
    cat(sprintf("  * Number of sites: %s\n", format(as.integer(length(x$pos)),
                                                    big.mark = ",")))
    invisible(NULL)
}


#' Print output from `vars_phylo`
#'
#'
#' @noRd
#' @export
#'
print.vars_phylo_info <- function(x, digits = max(3, getOption("digits") - 3), ...) {

    cat("< Phylo variant-creation info >\n")
    cat(sprintf("  * Number of variants: %i\n", length(x$phylo[[1]]$tip.label)))
    cat(sprintf("  * Number of trees: %i\n", length(x$phylo)))
    invisible(NULL)

}


#' Print output from `vars_theta`
#'
#'
#' @noRd
#' @export
#'
print.vars_theta_info <- function(x, digits = max(3, getOption("digits") - 3), ...) {

    cat("< Theta variant-creation info >\n")
    cat(sprintf("  * Theta: %.3g\n", x$theta))
    cat("# Phylogenetic tree:\n")
    print(x$phylo)
    invisible(NULL)

}


#' Print output from `vars_gtrees`
#'
#'
#'
#' @noRd
#' @export
#'
print.vars_gtrees_info <- function(x, digits = max(3, getOption("digits") - 3), ...) {

    cat("< Gene trees variant-creation info >\n")
    cat(sprintf("  * Number of chromosomes: %i\n", length(x$trees)))
    invisible(NULL)

}







# ====================================================================================`
# ====================================================================================`

# * TO_VAR_SET * -----

# ====================================================================================`
# ====================================================================================`

#' Used to convert info to VarSet pointer.
#'
#' @noRd
#'
to_var_set <- function (x, reference, sub, ins, del, epsilon, n_threads, show_progress) {
    UseMethod("to_var_set", x)
}

#' Create variants from segregating-site info from coalescent simulations.
#'
#'
#' @noRd
#'
to_var_set.vars_ssites_info <- function(x, reference, sub, ins, del, epsilon,
                                        n_threads, show_progress) {


    chrom_sizes <- reference$sizes()

    # Fill and check the position column in `x$mats`
    x$mats <- fill_coal_mat_pos(x$mats, chrom_sizes)

    variants_ptr <- add_ssites_cpp(reference$ptr(),
                                   x$mats,
                                   sub$Q(),
                                   sub$pi_tcag(),
                                   ins$rates(),
                                   del$rates(),
                                   n_threads,
                                   show_progress)

    return(variants_ptr)

}


#' Create variants from VCF file
#'
#'
#' @noRd
#'
to_var_set.vars_vcf_info <- function(x, reference, sub, ins, del, epsilon,
                                     n_threads, show_progress) {

    chrom_names <- view_ref_genome_chrom_names(reference$ptr())
    unq_chrom <- unique(x$chrom)

    if (!all(unq_chrom %in% chrom_names)) {
        if (x$print_names) print(unq_chrom)
        stop("\nChromosome name(s) in VCF file don't match those in the ",
             "`ref_genome` object. ",
             "It's probably easiest to manually change the `ref_genome` object ",
             "(using `$set_names()` method) to have the same names as the VCF file. ",
             "Re-run `vars_vcf` with `print_names = TRUE` to see the VCF-file names.",
             call. = FALSE)
    }

    # Converts items in `chrom` to 0-based indices of chromosomes in ref. genome
    chrom_inds <- match(x$chrom, chrom_names) - 1


    variants_ptr <- read_vcfr(reference$ptr(), x$var_names,
                              x$haps, chrom_inds, x$pos, x$ref_chrom)

    return(variants_ptr)

}


#' Create variants from phylogenetic tree(s).
#'
#'
#' @noRd
#'
to_var_set.vars_phylo_info <- function(x, reference, sub, ins, del, epsilon,
                                       n_threads, show_progress) {

    phy <- x$phylo
    class(phy) <- "list"

    n_vars <- length(phy$tip.label)
    n_chroms <- as.integer(reference$n_chroms())

    if (!length(phy) %in% c(1L, n_chroms)) {
        stop("\nIn function `vars_phylo`, you must provide information for 1 tree ",
             "or a tree for each reference genome chromosome. ",
             "It appears you need to re-run `vars_phylo` before attempting to ",
             "run `create_variants` again.")
    }

    if (length(phy) == 1) phy <- rep(phy, n_chroms)

    trees_ptr <- phylo_to_ptr(phy, n_chroms)

    var_set_ptr <- trees_to_var_set(trees_ptr, reference, sub, ins, del, epsilon,
                                    n_threads, show_progress)

    return(var_set_ptr)

}


#' Create variants from theta parameter.
#'
#'
#' @noRd
#'
to_var_set.vars_theta_info <- function(x,
                                       reference,
                                       sub, ins, del, epsilon,
                                       n_threads, show_progress) {

    phy <- x$phylo
    theta <- x$theta

    n_vars <- length(phy$tip.label)
    n_chroms <- reference$n_chroms()

    # Calculating L from theta:
    # E(L) = 4 * N * a; a = sum(1 / (1:(n_chroms-1)))
    a <- sum(1 / (1:(n_vars-1)))
    # theta = 4 * N * mu
    # ------------*
    # Calculating mu:
    # ------------*
    # Indel rates (same for each nucleotide):
    indel <- sum(ins$rates() * 0.25) + sum(del$rates() * 0.25)
    # Average mutation rate among all nucleotides:
    mu <- sum({rowSums(sub$Q()) + indel} * sub$pi_tcag())
    # ------------*
    # So if we know theta and mu, then...
    # ------------*
    L <- theta * a / mu
    # Now rescale to have total tree length of `L`:
    phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * L

    trees_ptr <- phylo_to_ptr(phy, n_chroms)

    var_set_ptr <- trees_to_var_set(trees_ptr, reference, sub, ins, del, epsilon,
                                    n_threads, show_progress)

    return(var_set_ptr)

}


#' Create variants from gene trees.
#'
#'
#' @noRd
#'
to_var_set.vars_gtrees_info <- function(x, reference, sub, ins, del, epsilon,
                                        n_threads, show_progress) {

    trees_ptr <- gtrees_to_ptr(x$trees, reference)

    var_set_ptr <- trees_to_var_set(trees_ptr, reference, sub, ins, del, epsilon,
                                    n_threads, show_progress)

    return(var_set_ptr)

}
