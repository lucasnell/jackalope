


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






