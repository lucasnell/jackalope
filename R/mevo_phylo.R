


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





#' Read info from a `phylo` object.
#'
#' @inheritParams make_phylo_info
#'
#' @return An external pointer to the phylogenetic info needed to do the sequence
#'     simulations.
#'
#' @noRd
#'
read_phy_obj <- function(phy, n_seqs, chunked, err_msg = "") {

    if (!inherits(phy, "phylo") & !inherits(phy, "multiPhylo") &
        !inherits(phy, "list")) {
        stop(sprintf(err_msg, "of class phylo, multiPhylo, or list"), call. = FALSE)
    }
    if (inherits(phy, "list")) {
        if (!all(sapply(phy, inherits, what = "phylo"))) {
            stop(paste0("\nFor `create_variants` method \"phy\", ",
                        "if providing a list, the list must only contain `phylo` ",
                        "objects."),
                 call. = FALSE)
        }
    }
    if (inherits(phy, "multiPhylo") | inherits(phy, "list")) {
        if (length(phy) != n_seqs) {
            stop(paste0("\nFor `create_variants` method \"phy\", ",
                        "if providing a list or `multiPhylo` object, its length ",
                        "must equal the number of sequences."),
                 call. = FALSE)
        }
    }
    # So all inputs are lists of the proper length:
    if (inherits(phy, "phylo")) phy <- rep(list(phy), n_seqs)
    if (inherits(phy, "multiPhylo")) class(phy) <- "list"

    phylo_info <- lapply(phy,
                         function(p) {
                             p <- ape::reorder.phylo(p, order = "cladewise")
                             labels <- paste(p$tip.label)  # <-- making sure they're strings
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







#' Read info from a coalescent object from scrm or coala.
#'
#' @inheritParams make_phylo_info
#'
#' @return An XPtr to the info needed from the phylogenies to do the sequence
#'     simulations.
#'
#' @noRd
#'
#'
read_coal_obj <- function(coal_obj, seq_sizes, chunked, err_msg) {

    # Check for coal_obj begin a list and either having a `trees` field or all its
    # items within having `trees` fields
    err <- FALSE
    nested <- FALSE
    if (!inherits(coal_obj, "list")) {
        err <- TRUE
    } else if (is.null(coal_obj$trees)) {
        if (any(sapply(coal_obj, function(x) is.null(x$trees)))) {
            err <- TRUE
        } else {
            nested <- TRUE
            if (any(sapply(coal_obj, function(x) length(x$trees)) != 1)) err <- TRUE
            if (any(!sapply(coal_obj, function(x) inherits(x$trees, "list")))) err <- TRUE
        }
    }
    if (err) {
        stop(sprintf(err_msg,
                     paste("(1) a list with a `trees` field present or",
                           "(2) a list of lists, each sub-list containing a `trees`",
                           "field of length 1. For more, see `?create_variants`.")),
             call. = FALSE)
    }

    if (nested) {
        trees <- lapply(coal_obj, function(x) x$trees[[1]])
    } else trees <- coal_obj$trees

    if (length(trees) != length(seq_sizes)) {
        stop(sprintf(err_msg,
                     paste("result in a number of trees that's the same as the number",
                           "of sequences.",
                           "For more, see `?create_variants`.")),
             call. = FALSE)
    }

    phylo_info <- mapply(process_coal_tree_string, trees, seq_sizes,
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)

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
read_ms_output <- function(ms_filename, seq_sizes, chunked, err_msg) {

    if (!single_string(ms_filename)) {
        stop(sprintf(err_msg, "a single string"), call. = FALSE)
    }

    trees <- read_ms_output_(ms_filename)

    phylo_info <- mapply(process_coal_tree_string, trees, seq_sizes,
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)

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




#' Read info from a NEWICK file.
#'
#' @inheritParams make_phylo_info
#'
#' @return An external pointer to the phylogenetic info needed to do the sequence
#'     simulations.
#'
#' @noRd
#'
read_newick <- function(newick_filename, n_seqs, chunked, err_msg) {

    if (!vec_string(newick_filename, c(1, n_seqs))) {
        stop(sprintf(err_msg, paste("a single string or a vector of strings of the",
                                    "same length as the number of sequences")),
             call. = FALSE)
    }

    phy <- lapply(newick_filename, ape::read.tree)

    if (length(phy) == 1) phy <- rep(phy, n_seqs)

    trees_ptr <- read_phy_obj(phy, n_seqs, chunked, err_msg)

    return(trees_ptr)
}



#' Random phylogenetic tree from theta and mu parameters.
#'
#' Note that mu should be derived from the mutation object, not passed to the function
#' by the user.
#'
#' @inheritParams make_phylo_info
#'
#' @return An external pointer to the phylogenetic info needed to do the sequence
#'     simulations.
#'
#' @noRd
#'
read_theta <- function(theta_n_vars, mu, n_seqs, chunked, err_msg) {

    err_msg <- sprintf(err_msg,
                       paste("a named list or numeric vector, with the names",
                             "\"theta\" and \"n_vars\".",
                             "\"theta\" must be single number, and",
                             "\"n_vars\" must be a single whole number."))
    if ((!inherits(theta_n_vars, "list") & !inherits(theta_n_vars, "numeric")) |
        is.null(names(theta_n_vars))) {
        stop(err_msg, call. = FALSE)
    }
    if (!single_number(theta_n_vars[["theta"]])) stop(err_msg, call. = FALSE)
    if (!single_whole_number(theta_n_vars[["n_vars"]])) stop(err_msg, call. = FALSE)

    theta <- theta_n_vars[["theta"]]
    n_vars <- theta_n_vars[["n_vars"]]

    # Generate random coalescent tree:
    phy <- ape::rcoal(n_vars)

    # Calculating L from theta:
    # E(L) = 4 * N * a; a = sum(1 / (1:(n_seqs-1)))
    a <- sum(1 / (1:(n_vars-1)))
    # theta = 4 * N * mu
    # So if we know theta and mu, then...
    L <- theta * a / mu

    # Now rescale to have total tree length of `L`:
    phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * L

    trees_ptr <- read_phy_obj(phy, n_seqs, chunked)

    return(trees_ptr)

}




#' Create phylogenetic information object from one of multiple methods.
#'
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
make_phylo_info <- function(method,
                            method_info,
                            seq_sizes,
                            n_seqs,
                            mu,
                            chunk_size) {

    chunked <- chunk_size > 0

    err_msg <- paste0("\nFor `create_variants` method \"", method,
                      "\", `method_info` must be %s.")


    if (method == "phylo") {

        trees_ptr <- read_phy_obj(method_info, n_seqs, chunked, err_msg)

    } else if (method == "coal_obj") {

        trees_ptr <- read_coal_obj(method_info, seq_sizes, chunked, err_msg)

    } else if (method == "ms_file") {

        trees_ptr <- read_ms_output(method_info, seq_sizes, chunked, err_msg)

    } else if (method == "newick") {

        trees_ptr <- read_newick(method_info, n_seqs, chunked, err_msg)

    } else if (method == "theta") {

        trees_ptr <- read_theta(method_info, mu, n_seqs, chunked, err_msg)

    } else {
        stop("\nAn improper method argument was input to the `make_phylo_info` function",
             call. = FALSE)
    }

    return(trees_ptr)
}






