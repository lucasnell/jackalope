
#' Process one gene-tree string from a coalescent simulator with ms-style output.
#'
#' @param str The string to process.
#'
#' @noRd
#'
process_coal_tree_string <- function(str) {

    if (length(str) > 1) {
        if (!all(grepl("^\\[", str))) {
            stop("\nThe coalescent string appears to include ",
                 "recombination but does not include sizes for each region.",
                 call. = FALSE)
        }
        sizes_ <- as.numeric(sapply(str, function(x) strsplit(x, "\\[|\\]")[[1]][2]))
    } else {
        sizes_ <- 1
    }
    ends <- cumsum(sizes_) - 1
    starts <- c(0, head(ends, -1) + 1)

    phylo_ <- ape::read.tree(text = str)
    # If no recombination (so only one phylo per sequence),
    # then we need to make sure that it has the same nestedness as if there were
    # >1 phylo objects:
    if (inherits(phylo_, "phylo")) {
        phylo_ <- list(phylo_)
    }

    out <- rep(list(NA), length(phylo_))

    for (i in 1:length(phylo_)) {
        phy <- reorder(phylo_[[i]], order = "cladewise")
        labels <- paste(phy$tip.label) # used paste to make sure they're characters
        branch_lens <- phy$edge.length
        edges <- phy$edge
        out[[i]] <- list(branch_lens = branch_lens, edges = edges, labels = labels,
                         start = starts[i], end = ends[i])
    }

    return(out)
}






#' Read info from a coalescent object from scrm or coala.
#'
#' @param coal_obj The coalescent-simulation object.
#'
#' @return An XPtr to the info needed from the phylogenies to do the sequence simulations.
#'
#' @noRd
#'
#'
read_coal_obj <- function(coal_obj) {

    if (is.null(coal_obj$trees) | !inherits(coal_obj, "list")) {
        stop("\nWhen reading trees from a coalescent object from the scrm or coala ",
             "packages, the object must be a list with a `trees` field present. ",
             "In coala, make sure to add `+ sumstat_trees()` to the coalmodel, ",
             "and in scrm, make sure to add the `-T` option.",
             call. = FALSE)
    }

    trees <- coal_obj$trees

    tree_info <- lapply(trees, process_coal_tree_string)

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

    trees_ptr <- phylo_info_to_trees(tree_info)

    return(trees_ptr)
}




#' Read info from a NEWICK file.
#'
#' @param newick_filename The filename for the NEWICK phylogeny.
#' @param n_seqs The number of sequences in the reference genome.
#'
#' @return An XPtr to the info needed from the phylogenies to do the sequence simulations.
#'
#' @noRd
#'
read_newick <- function(newick_filename, n_seqs) {

    phy <- ape::read.tree(file = newick_filename)

    phy <- reorder(phy, order = "cladewise")
    labels <- paste(phy$tip.label) # used paste to make sure they're characters
    branch_lens <- phy$edge.length
    edges <- phy$edge
    phy_info <- list(branch_lens = branch_lens, edges = edges, labels = labels,
                     start = 0, end = 0)

    # I'm repeating this information to make one set of phylo. info per sequence
    # I'm doing `list(list(...))` to keep nestedness the same among the different methods
    tree_info <- rep(list(list(phy_info)), n_seqs)

    trees_ptr <- phylo_info_to_trees(tree_info)

    return(trees_ptr)
}



read_ms_output <- function(ms_filename) {

}
