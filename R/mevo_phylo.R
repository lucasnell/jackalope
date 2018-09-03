


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
            stop("\nThe coalescent string appears to include ",
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
#' @param phy The `phylo` object.
#' @param n_seqs The number of sequences in the reference genome.
#' @param chunked Boolean for whether the sampling for mutation locations
#'     will be done in chunks.
#'
#' @return An XPtr to the info needed from the phylogenies to do the sequence simulations.
#'
#' @noRd
#'
read_phy_obj <- function(phy, n_seqs, chunked) {

    phy <- ape::reorder.phylo(phy, order = "cladewise")
    labels <- paste(phy$tip.label) # used paste to make sure they're characters
    branch_lens <- phy$edge.length
    edges <- phy$edge
    phy_info <- list(branch_lens = branch_lens, edges = edges, labels = labels,
                     start = 0, end = 0)

    # I'm repeating this information to make one set of phylo. info per sequence
    # I'm doing `list(list(...))` to keep nestedness the same among the different methods
    tree_info <- rep(list(list(phy_info)), n_seqs)

    if (!chunked) {
        trees_ptr <- phylo_info_to_trees(tree_info)
    } else {
        trees_ptr <- phylo_info_to_trees_chunk(tree_info)
    }

    return(trees_ptr)
}







#' Read info from a coalescent object from scrm or coala.
#'
#' @param coal_obj The coalescent-simulation object.
#' @param seq_sizes Vector of sequence sizes.
#' @inheritParams read_phy_obj
#'
#' @return An XPtr to the info needed from the phylogenies to do the sequence simulations.
#'
#' @noRd
#'
#'
read_coal_obj <- function(coal_obj, seq_sizes, chunked) {

    if (is.null(coal_obj$trees) | !inherits(coal_obj, "list")) {
        stop("\nWhen reading trees from a coalescent object from the scrm or coala ",
             "packages, the object must be a list with a `trees` field present. ",
             "In coala, make sure to add `+ sumstat_trees()` to the coalmodel, ",
             "and in scrm, make sure to add the `-T` option.",
             call. = FALSE)
    }

    trees <- coal_obj$trees

    tree_info <- mapply(process_coal_tree_string, trees, seq_sizes,
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)

    # Making sure all labels are the same
    label_mat <- do.call(rbind,
                         lapply(tree_info, function(x) {
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
        trees_ptr <- phylo_info_to_trees(tree_info)
    } else {
        trees_ptr <- phylo_info_to_trees_chunk(tree_info)
    }

    return(trees_ptr)
}





#' Read info from ms-style output file.
#'
#' @param ms_filename The filename for ms-style output.
#' @inheritParams read_coal_obj
#' @inheritParams read_phy_obj
#'
#' @return An XPtr to the info needed from the phylogenies to do the sequence simulations.
#'
#' @noRd
#'
#'
read_ms_output <- function(ms_filename, seq_sizes, chunked) {

    trees <- read_ms_output_(ms_filename)

    tree_info <- mapply(process_coal_tree_string, trees, seq_sizes,
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)

    # Making sure all labels are the same
    label_mat <- do.call(rbind,
                         lapply(tree_info, function(x) {
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
        trees_ptr <- phylo_info_to_trees(tree_info)
    } else {
        trees_ptr <- phylo_info_to_trees_chunk(tree_info)
    }

    return(trees_ptr)

}




#' Read info from a NEWICK file.
#'
#' @param newick_filename The filename for the NEWICK phylogeny.
#' @inheritParams read_phy_obj
#'
#' @return An XPtr to the info needed from the phylogenies to do the sequence simulations.
#'
#' @noRd
#'
read_newick <- function(newick_filename, n_seqs, chunked) {

    phy <- ape::read.tree(file = newick_filename)

    trees_ptr <- read_phy_obj(phy, n_seqs, chunked)

    return(trees_ptr)
}



#' Random phylogenetic tree from theta and mu parameters.
#'
#' Note that mu should be derived from the mutation object, not passed to the function
#' by the user.
#'
#' @param theta Theta parameter, population-scaled mutation rate.
#' @param mu Average mutation rate (per bp per generation).
#' @inheritParams read_phy_obj
#'
#' @return An XPtr to the info needed from the phylogenies to do the sequence simulations.
#'
#' @noRd
#'
read_theta <- function(theta, mu, n_vars, n_seqs, chunked) {

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

