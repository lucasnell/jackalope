



# ======================================================================================`
# ======================================================================================`

# HELPER FUNCTIONS -----

# ======================================================================================`
# ======================================================================================`


#' Used to convert info to pointer to VarSet info.
#'
#' @noRd
#'
to_var_set <- function (x, reference, mevo_obj, n_threads, show_progress) {
    UseMethod("to_var_set", x)
}
#' Go from pointer to trees info to a pointer to a VarSet object
#'
#' Used below in `to_var_set` methods
#'
#' @noRd
#'
trees_to_var_set <- function(phylo_info_ptr, reference, mevo_obj, n_threads,
                             show_progress) {

    # Make sampler_base_ptr
    sampler_base_ptr <- mevo_obj_to_ptr(mevo_obj)

    # Make Gamma matrices (for mutation-rate variability among sites):
    gamma_mats <- mevo_obj$gamma_mats

    # Make variants pointer:
    variants_ptr <- NULL
    if (mevo_obj$chunk_size > 0) {
        variants_ptr <- evolve_seqs_chunk(
            reference$genome,
            sampler_base_ptr,
            phylo_info_ptr,
            gamma_mats,
            n_threads,
            show_progress)
    } else {
        variants_ptr <- evolve_seqs(
            reference$genome,
            sampler_base_ptr,
            phylo_info_ptr,
            gamma_mats,
            n_threads,
            show_progress)
    }

    return(variants_ptr)

}




#' Read info from a `phylo` object.
#'
#' @inheritParams make_phylo_info
#'
#' @return An external pointer to the phylogenetic info needed to do the sequence
#'     simulations.
#'
#' @noRd
#'
phylo_to_ptr <- function(phy, n_seqs, chunked) {

    if ((!inherits(phy, "phylo") && !inherits(phy, "multiPhylo") &&
        !inherits(phy, "list")) ||
        (inherits(phy, "list") && !all(sapply(phy, inherits, what = "phylo")))) {
        stop("\nThe `phy` argument to the internal function `phylo_to_ptr` should ",
             "only ever be an object of class \"phylo\", \"multiPhylo\", or a ",
             "list of \"phylo\" objects.")
    }

    if ((inherits(phy, "multiPhylo") || inherits(phy, "list")) &&
        length(phy) != n_seqs) {
        stop("\nThe `phy` argument to the internal function `phylo_to_ptr` should ",
             "have a length of 1 or equal to the number of sequences if it's of class",
                "\"multiPhylo\" or list.")
    }
    # So all inputs are lists of the proper length:
    if (inherits(phy, "phylo")) phy <- rep(list(phy), n_seqs)
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

    if (!chunked) {
        trees_ptr <- phylo_info_to_trees(phylo_info)
    } else {
        trees_ptr <- phylo_info_to_trees_chunk(phylo_info)
    }

    return(trees_ptr)
}




# ======================================================================================`
# ======================================================================================`

# GENE_TREES METHOD -----

# ======================================================================================`
# ======================================================================================`


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
#'
#' If using an output file from a command-line program like `ms`/`msms`,
#' add the `-T` option.
#'
#'
#' @param obj Object containing gene trees.
#'     This can be one of the following:
#'     (1) A single `list` with a `trees` field inside. This field must
#'     contain a set of gene trees for each sequence.
#'     (2) A list of lists, each sub-list containing a `trees` field of
#'     length 1. The top-level list must be of the same length as the
#'     number of sequences.
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
    type <- NULL

    if (!is.null(obj)) {

        type <- "object"

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

        type <- "file"

        if (!is_type(fn, "character", 1)) {
            err_msg("vars_gtrees", "fn", "NULL or a single string")
        }

        trees <- read_ms_trees_(fn)

        if (any(sapply(trees, length) == 0)) {
            stop("\nIn ms-style output file, one or more sequences have no trees.",
                 call. = FALSE)
        }

    }

    out <- list(trees = trees, type = type)
    class(out) <- "vars_gtrees_info"

    return(out)

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
    cat(sprintf("  * Number of sequences: %i\n", length(x$trees)))
    invisible(NULL)

}


to_var_set.vars_gtrees_info <- function(x, reference, mevo_obj,
                                        n_threads, show_progress) {

    trees_ptr <- NULL

    if (x$type == "object") {
        trees_ptr <- read_coal_trees(x$trees, reference, mevo_obj)
    } else if (x$type == "file") {
        trees_ptr <- read_ms_trees(x$trees, reference, mevo_obj)
    } else {
        stop("\nAll `vars_gtrees_info` objects should have a type field that's either ",
             "\"object\" or \"file\".", call. = FALSE)
    }

    var_set_ptr <- trees_to_var_set(trees_ptr, reference, mevo_obj, n_threads,
                                    show_progress)

    return(var_set_ptr)

}


#' Process one gene-tree string from a coalescent simulator with ms-style output.
#'
#' @param str The string to process.
#' @param seq_size The number of bp in the sequence associated with the input string.
#'
#' @noRd
#'
process_coal_tree_string <- function(str, seq_size) {

    if (length(str) > 1) {
        if (!all(grepl("^\\[", str))) {
            stop("\nA coalescent string appears to include ",
                 "recombination but does not include sizes for each region.",
                 call. = FALSE)
        }
        sizes_ <- as.numeric(sapply(str, function(x) strsplit(x, "\\[|\\]")[[1]][2]))
        # If they're <= 1, then they're not # bp, they're proportion of sequence
        if (all(sizes_ <= 1) & seq_size > 1) {
            sizes_ <- sizes_ / sum(sizes_)
            sizes_ <- round(sizes_ * seq_size, 0)
            # Remove any zero sizes:
            sizes_ <- sizes_[sizes_ > 0]
            # If there's nothing left, just make it of length 1:
            if (length(sizes_) == 0) {
                sizes_ <- seq_size
            # If it doesn't round quite right, then randomly add/subtract:
            } else if (sum(sizes_) != seq_size) {
                inds <- sample.int(length(sizes_), abs(seq_size - sum(sizes_)))
                sizes_[inds] <- sizes_[inds] + sign(seq_size - sum(sizes_))
            }
        } else if (sum(sizes_) != seq_size) {
            stop("\nA coalescent string appears to include ",
                 "recombination but the combined sizes of all regions don't match ",
                 "the size of the sequence.", call. = FALSE)
        }
    } else {
        sizes_ <- 1
    }
    ends <- cumsum(sizes_) - 1
    starts <- c(0, utils::head(ends, -1) + 1)

    phylo_ <- ape::read.tree(text = str)
    # If no recombination (so only one phylo per sequence),
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

#' Read gene-tree info from a coalescent object from scrm or coala.
#'
#' @inheritParams make_phylo_info
#'
#' @return An XPtr to the info needed from the phylogenies to do the sequence
#'     simulations.
#'
#' @noRd
#'
#'
read_coal_trees <- function(trees, reference, mevo_obj) {

    seq_sizes <- reference$sizes()
    chunked <- mevo_obj$chunk_size > 0

    if (length(trees) != length(seq_sizes)) {
        err_msg("create_variants", "method_info",
                "a list that results in a number of trees that's the same as the",
                "number of sequences in the reference genome,",
                "when `method` = \"coal_trees\"")
    }

    phylo_info <- mapply(process_coal_tree_string, trees, seq_sizes,
                         SIMPLIFY = FALSE, USE.NAMES = FALSE)

    unq_n_tips <- lapply(phylo_info,
                         function(x) sapply(x, function(xx) length(xx$labels)))
    unq_n_tips <- unique(do.call(c, unq_n_tips))
    if (length(unq_n_tips) > 1) {
        stop("\nIn the input coalescent object, the gene trees don't all have the ",
             "same number of tips.", call. = FALSE)
    }
    unq_tips_names <- sapply(phylo_info,
                             function(x) {
                                 tips_ <- do.call(c, lapply(x, function(xx) xx$labels))
                                 paste(sort(unique(tips_)), collapse = "___")
                             })
    if (length(unique(unq_tips_names)) > 1) {
        stop("\nIn the input coalescent file, the gene trees don't all have the ",
             "same tip names.", call. = FALSE)
    }

    # Making sure all labels are the same
    label_mat <- do.call(rbind,
                         lapply(phylo_info, function(x) {
                             t(sapply(x, function(xx) xx$labels))
                         }))
    label_mat <- t(apply(label_mat, 1, sort))
    for (i in 2:nrow(label_mat)) {
        if (any(label_mat[1,] != label_mat[i,])) {
            stop("\nIn the input coalescent object, not all labels are the same.",
                 call. = FALSE)
        }
    }

    if (!chunked) {
        trees_ptr <- phylo_info_to_trees(phylo_info)
    } else {
        trees_ptr <- phylo_info_to_trees_chunk(phylo_info)
    }

    return(trees_ptr)
}







#' Read info from ms-style output file.
#'
#' @inheritParams make_phylo_info
#'
#' @return An external pointer to the phylogenetic info needed to do the sequence
#'     simulations.
#'
#' @noRd
#'
#'
read_ms_trees <- function(trees, reference, mevo_obj) {

    seq_sizes <- reference$sizes()
    chunked <- mevo_obj$chunk_size > 0

    phylo_info <- mapply(process_coal_tree_string, trees, seq_sizes,
                         SIMPLIFY = FALSE, USE.NAMES = FALSE)

    unq_n_tips <- lapply(phylo_info,
                         function(x) sapply(x, function(xx) length(xx$labels)))
    unq_n_tips <- unique(do.call(c, unq_n_tips))
    if (length(unq_n_tips) > 1) {
        stop("\nIn the input coalescent file, the gene trees don't all have the ",
             "same number of tips.", call. = FALSE)
    }
    unq_tips_names <- sapply(phylo_info,
                             function(x) {
                                 tips_ <- do.call(c, lapply(x, function(xx) xx$labels))
                                 paste(sort(unique(tips_)), collapse = "___")
                             })
    if (length(unique(unq_tips_names)) > 1) {
        stop("\nIn the input coalescent file, the gene trees don't all have the ",
             "same tip names.", call. = FALSE)
    }


    # Making sure all labels are the same
    label_mat <- do.call(rbind,
                         lapply(phylo_info, function(x) {
                             t(sapply(x, function(xx) xx$labels))
                         }))
    label_mat <- t(apply(label_mat, 1, sort))
    for (i in 2:nrow(label_mat)) {
        if (any(label_mat[1,] != label_mat[i,])) {
            stop("\nIn the input coalescent object, not all labels are the same.",
                 call. = FALSE)
        }
    }

    if (!chunked) {
        trees_ptr <- phylo_info_to_trees(phylo_info)
    } else {
        trees_ptr <- phylo_info_to_trees_chunk(phylo_info)
    }

    return(trees_ptr)

}




# ======================================================================================`
# ======================================================================================`

# PHYLO METHOD -----

# ======================================================================================`
# ======================================================================================`


#' Create necessary information to create variants using phylogenetic tree(s)
#'
#' This function organizes higher-level information for creating variants from
#' phylogenetic tree(s) output as `phylo` objects or NEWICK files.
#'
#' @param obj Object containing phylogenetic tree(s).
#'     This can be (1) a single \code{\link[ape]{phylo}} object that represents all
#'     sequences in the genome or
#'     (2) a `list` or `multiPhylo` object containing a `phylo` object for
#'     each reference sequence.
#'     In the latter case, phylogenies will be assigned to sequences in the
#'     order provided.
#'     Defaults to `NULL`.
#' @param fn One or more string(s), each of which specifies the file name
#'     of a NEWICK file containing a phylogeny.
#'     If one name is provided, that phylogeny will be used for all sequences.
#'     If more than one is provided, there must be a phylogeny for each reference
#'     genome sequence, and phylogenies will be assigned to sequences
#'     in the order provided.
#'     Defaults to `NULL`.
#'
#'
#' @return A `vars_phylo_info` object containing information used in `create_variants`
#'     to create haploid variants.
#'     This class is just a wrapper around a list containing phylogenetic tree
#'     information for each reference sequence.
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
                    "of class \"phylo\", \"multiPhylo\", or a list of \"phylo\" objects")
        }

        phy <- obj
        # So all output is a list
        if (inherits(phy, "phylo")) phy <- list(phy)
        if (inherits(phy, "multiPhylo")) class(phy) <- "list"

    }
    if (!is.null(fn)) {
        if (!is_type(fn, "character")) {
            err_msg("vars_phylo", "fn", "a character vector")
        }
        phy <- lapply(fn, ape::read.tree)
    }

    out <- list(phylo = phy)
    class(out) <- "vars_phylo_info"

    return(out)

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

to_var_set.vars_phylo_info <- function(x, reference, mevo_obj,
                                       n_threads, show_progress) {

    n_vars <- length(phy$tip.label)
    n_seqs <- as.integer(reference$n_seqs())
    chunked <- mevo_obj$chunk_size > 0

    phy <- x$phylo
    class(phy) <- "list"

    if (!length(phy) %in% c(1L, n_seqs)) {
        stop("\nIn function `vars_phylo`, you must provide information for 1 tree ",
             "or a tree for each reference genome sequence. ",
             "It appears you need to re-run `vars_phylo` before attempting to ",
             "run `create_variants` again.")
    }

    if (length(phy) == 1) phy <- rep(phy, n_seqs)

    trees_ptr <- phylo_to_ptr(phy, n_seqs, chunked)

    var_set_ptr <- trees_to_var_set(trees_ptr, reference, mevo_obj, n_threads,
                                    show_progress)

    return(var_set_ptr)

}




# ======================================================================================`
# ======================================================================================`

# THETA METHOD -----

# ======================================================================================`
# ======================================================================================`


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

to_var_set.vars_theta_info <- function(x, reference, mevo_obj,
                                       n_threads, show_progress) {

    phy <- x$phylo
    theta <- x$theta

    n_vars <- length(phy$tip.label)
    n_seqs <- reference$n_seqs()
    chunked <- mevo_obj$chunk_size > 0

    # Calculating L from theta:
    # E(L) = 4 * N * a; a = sum(1 / (1:(n_seqs-1)))
    a <- sum(1 / (1:(n_vars-1)))
    # theta = 4 * N * mu
    # So if we know theta and mu, then...
    L <- theta * a / mevo_obj$mu()
    # Now rescale to have total tree length of `L`:
    phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * L

    trees_ptr <- phylo_to_ptr(phy, n_seqs, chunked)

    var_set_ptr <- trees_to_var_set(trees_ptr, reference, mevo_obj, n_threads,
                                    show_progress)

    return(var_set_ptr)

}








# ======================================================================================`
# ======================================================================================`

# COMBINED -----

# ======================================================================================`
# ======================================================================================`

# This function will probably need to be removed or changed drastically


#' Create phylogenetic information object from one of multiple methods.
#'
#'
#' \describe{
#'     \item{`method = "phylo"`}{One of the following object types is allowed:
#'         \itemize{
#'             \item A single \code{\link[ape]{phylo}} object that represents all
#'                 sequences in the genome.
#'             \item A `list` or `multiPhylo` object containing a `phylo` object for
#'                 each reference sequence.
#'                 Phylogenies will be assigned to sequences in the order provided.
#'             \item One or more string(s), each of which specifies
#'                 a name of a NEWICK file containing a phylogeny.
#'                 If one name is provided, that phylogeny will be used for
#'                 all sequences.
#'                 If more than one is provided, there must be a phylogeny for
#'                 each sequence, and phylogenies will be assigned to sequences
#'                 in the order provided.
#'         }
#'     }
#'     \item{`method = "coal_trees"`}{One of the following object types is allowed:
#'         \itemize{
#'             \item A single `list` with a `trees` field inside. This field must
#'                 contain a set of gene trees for each sequence.
#'             \item A list of lists, each sub-list containing a `trees` field of
#'                 length 1. The top-level list must be of the same length as the
#'                 number of sequences.
#'             \item A single string specifying the name of the file containing
#'                 the `ms`-style coalescent output with gene trees.
#'         }
#'         The top two options are designed after the `trees` fields in the output from
#'         the `scrm` and `coala` packages.
#'         (These packages are not required to be installed when installing
#'         `jackalope`.)
#'         To get gene trees, make sure to add `+ sumstat_trees()`
#'         to the `coalmodel` for `coala`, or
#'         make sure that `"-T"` is present in `args` for `scrm`.
#'         If using an output file from a command-line program like `ms`/`msms`,
#'         add the `-T` option.
#'     }
#'     \item{`method = "theta"`}{A named vector or list containing the fields `theta`
#'         and `n_vars`, specifying the theta parameter (population-scaled mutation rate)
#'         and number of desired variants, respectively.}
#' }
#'
#' @inheritParams create_variants
#' @param seq_sizes Vector of sequence sizes.
#' @param n_seqs The number of sequences in the reference genome.
#' @param mu Average mutation rate (per bp per generation).
#'
#' @return An external pointer to the phylogenetic info needed to do the sequence
#'     simulations.
#'
#' @noRd
#'
make_phylo_info <- function(obj,
                            reference,
                            mevo_obj) {

    chunk_size <- mevo_obj$chunk_size
    seq_sizes <- reference$sizes()
    mu <- mevo_obj$mu()
    n_seqs <- length(seq_sizes)

    chunked <- chunk_size > 0

    if (method == "phylo") {

        trees_ptr <- phylo_to_ptr(method_info, n_seqs, chunked)

    } else if (method == "coal_trees") {

        if (is_type(method_info, "character", 1)) {
            trees_ptr <- read_ms_trees(method_info, seq_sizes, chunked)
        } else if (inherits(method_info, "list")) {
            trees_ptr <- read_coal_trees(method_info, seq_sizes, chunked)
        } else {
            err_msg("create_variants", "method_info", "a single string or a list",
                    "when `method` = \"coal_trees\"")
        }

    } else if (method == "newick") {

        trees_ptr <- read_newick(method_info, n_seqs, chunked)

    } else if (method == "theta") {

        trees_ptr <- read_theta(method_info, mu, n_seqs, chunked)

    } else {
        err_msg("make_phylo_info", "method", "one of \"phylo\", \"coal_trees\",",
                "\"newick\", or \"theta\" for a phylogenetic method")
    }

    return(trees_ptr)
}






