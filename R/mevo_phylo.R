


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


