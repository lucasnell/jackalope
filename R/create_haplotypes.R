
# ====================================================================================`
# ====================================================================================`

# * HELPERS * -----

# ====================================================================================`
# ====================================================================================`


# -------------*
#  Non-phylogenomic -----
# -------------*


#' Check validity of position columns in segregating-sites matrices.
#'
#' Used in `to_hap_set` for the ssites method.
#'
#' @noRd
#'
fill_coal_mat_pos <- function(sites_mats, chrom_sizes) {

    if (length(sites_mats) != length(chrom_sizes)) {
        stop("\nIn function `haps_ssites`, there must be exactly one segregating sites ",
             "matrix for each reference genome chromosome. ",
             "It appears you need to re-run `haps_ssites` before attempting to ",
             "run `create_haplotypes` again.")
    }

    for (i in 1:length(sites_mats)) {
        if (nrow(sites_mats[[i]]) == 0) next;
        pos <- sites_mats[[i]][,1]
        if (all(pos < 1 & pos > 0)) {
            # Converting to integer positions (0-based):
            pos  <- as.integer(pos * chrom_sizes[i]);
            sites_mats[[i]][,1] <- pos
            # There might be some repeats now, so removing those:
            if (anyDuplicated(pos) != 0) {
                dups <- duplicated(pos)
                sites_mats[[i]] <- sites_mats[[i]][!dups,]
            }
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
                     "It appears you need to re-run `haps_ssites` before attempting to ",
                     "run `create_haplotypes` again.")
            }
            sites_mats[[i]][,1] <- pos
        }

    }

    return(sites_mats)

}




# -------------*
#  Phylogenomic -----
# -------------*


# These first two are helpers for multiple phylogenomic methods


#' Go from pointer to trees info to a pointer to a VarSet object
#'
#' Used below in theta, phylo, and gtrees `to_hap_set` methods
#'
#' @noRd
#'
trees_to_hap_set <- function(trees_info, reference, sub, ins, del, epsilon,
                             n_threads, show_progress) {

    haplotypes_ptr <- evolve_across_trees(reference$ptr(),
                                        trees_info,
                                        sub$Q(),
                                        sub$U(),
                                        sub$Ui(),
                                        sub$L(),
                                        sub$invariant(),
                                        ins$rates(),
                                        del$rates(),
                                        epsilon,
                                        sub$pi_tcag(),
                                        n_threads,
                                        show_progress)

    return(haplotypes_ptr)

}





#' Process phylogenetic trees.
#'
#' It also standardizes tip indices so that in all phylogenies, edge-matrix indices refer
#' to the same tips.
#'
#' It does NOT create a sensible `n_bases` field!
#'
#' Used in `phylo_to_info_list` and `process_coal_tree_string`
#'
#' @noRd
#'
process_phy <- function(phy, ordered_tip_labels) {

    if (!ape::is.binary.phylo(phy)) {
        stop("\nAll phylogenetic trees must be binary. An option to remedy this might ",
             "be the function `ape::multi2di`.", call. = FALSE)
    }
    if (!ape::is.rooted(phy)) {
        stop("\nAll phylogenetic trees must be rooted. An option to remedy this might ",
             "be the function `ape::root`.", call. = FALSE)
    }

    # Order phylogeny so that extra objects at nodes can be cleared away ASAP:
    phy <- ape::reorder.phylo(phy, order = "cladewise")

    # Make sure tip labels are strings:
    phy$tip.label <- paste(phy$tip.label)

    # -----------*
    # Standarize tips:
    # -----------*

    if (!identical(sort(ordered_tip_labels), sort(phy$tip.label))) {
        stop("\nOne or more trees have differing tip labels.", call. = FALSE)
    }

    # If tips need re-ordering, do that now:
    if (!identical(ordered_tip_labels, phy$tip.label)) {

        new_phy <- phy

        m <- match(ordered_tip_labels, phy$tip.label)

        for (i in 1:length(m)) {
            new_phy$edge[,2][phy$edge[,2] == m[i]] <- i
        }

        new_phy$tip.label <- ordered_tip_labels

        phy <- new_phy
    }


    # -----------*
    # Use tips as nodes (so no intermediate objects must be made):
    # -----------*
    n_tips <- length(phy$tip.label)
    n_nodes <- phy$Nnode
    edge <- phy$edge

    for (i in 1:n_nodes) {
        edge_i <- edge[edge[,1] == i + n_tips, ]
        edge_i[,1] <- utils::tail(edge_i[,2], 1)
        edge[edge[,1] == i + n_tips, ] <- edge_i
        edge[edge[,2] == i + n_tips, 2] <- utils::tail(edge_i[,2], 1)
    }
    phy$edge <- edge


    # -----------*
    # Extract all info:
    # -----------*

    phy_info <- list(branch_lens = phy$edge.length,
                     edges = phy$edge,
                     labels = phy$tip.label,
                     n_bases = 0)

    return(phy_info)

}




#' Organize info into list to later create C++ tree class from phylo info.
#'
#' Used in `to_hap_set` for phylo and theta methods.
#'
#' @return An external pointer to the phylogenetic info needed to do the simulations.
#'
#' @noRd
#'
phylo_to_info_list <- function(phy, reference) {

    chrom_sizes <- reference$sizes()

    if (!inherits(phy, "list") || !all(sapply(phy, inherits, what = "phylo"))) {
        stop("\nThe `phy` argument to the internal function `phylo_to_info_list` should ",
             "only ever be a list of \"phylo\" objects.")
    }

    # Ordered tip labels:
    otl <- phy[[1]]$tip.label

    # Phylogeny information:
    phylo_info <- lapply(phy, process_phy, ordered_tip_labels = otl)
    for (i in 1:length(phylo_info)) phylo_info[[i]][["n_bases"]] <- chrom_sizes[i]

    # For proper nestedness:
    phylo_info <- lapply(phylo_info, function(x) list(x))

    return(phylo_info)
}



# __gtrees -----


#' Process one gene-tree string from a coalescent simulator with ms-style output.
#'
#' Used in `gtrees_to_info_list`.
#'
#' @param str The string to process.
#' @param chrom_size The number of bp in the chromosome associated with the input string.
#'
#' @noRd
#'
process_coal_tree_string <- function(str, chrom_size, ordered_tip_labels) {

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
        sizes_ <- chrom_size
    }


    phylo_ <- ape::read.tree(text = str)
    # If no recombination (so only one phylo per chromosome),
    # then we need to make sure that it has the same nestedness as if there were
    # >1 phylo objects:
    if (inherits(phylo_, "phylo")) {
        phylo_ <- list(phylo_)
    }

    out <- lapply(phylo_, process_phy, ordered_tip_labels = ordered_tip_labels)

    for (i in 1:length(phylo_)) out[[i]][["n_bases"]] <- sizes_[i]

    return(out)
}


#' Organize info into list to later create C++ tree class from gene-tree info.
#'
#' Used in `to_hap_set` for gtrees method.
#'
#' @return An XPtr to the info needed from the gene trees to do the simulations.
#'
#' @noRd
#'
#'
gtrees_to_info_list <- function(trees, reference) {

    chrom_sizes <- reference$sizes()

    if (length(trees) != length(chrom_sizes)) {
        stop("\nFor the gene-trees method of haplotype creation, there must be a set ",
             "of gene trees for each reference genome chromosome. ",
             "It appears you need to re-run `haps_gtrees` before attempting to ",
             "run `create_haplotypes` again.")
    }

    otl <- paste(ape::read.tree(text = trees[[1]][1])[["tip.label"]])

    # The process_phy function inside process_coal_tree_string does checking of tip names
    phylo_info <- mapply(process_coal_tree_string, trees, chrom_sizes,
                         MoreArgs = list(ordered_tip_labels = otl),
                         SIMPLIFY = FALSE, USE.NAMES = FALSE)

    return(phylo_info)
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
to_hap_set <- function(x, reference, sub, ins, del, epsilon, n_threads, show_progress) {

    fun <- NULL

    if (inherits(x, "haps_vcf_info")) {
        fun <- to_hap_set__haps_vcf_info
    } else if (inherits(x, "haps_ssites_info")) {
        fun <- to_hap_set__haps_ssites_info
    } else if (inherits(x, "haps_theta_info")) {
        fun <- to_hap_set__haps_theta_info
    } else if (inherits(x, "haps_phylo_info")) {
        fun <- to_hap_set__haps_phylo_info
    } else if (inherits(x, "haps_gtrees_info")) {
        fun <- to_hap_set__haps_gtrees_info
    } else stop("Unknown input to `haps_info` arg in `create_haplotypes`")

    haplotypes_ptr <- fun(x = x, reference = reference,
                        sub = sub, ins = ins, del = del, epsilon = epsilon,
                        n_threads = n_threads, show_progress = show_progress)

    return(haplotypes_ptr)

}

#' Create haplotypes from segregating-site info from coalescent simulations.
#'
#'
#' @noRd
#'
to_hap_set__haps_ssites_info <- function(x, reference, sub, ins, del, epsilon,
                                        n_threads, show_progress) {


    chrom_sizes <- reference$sizes()

    # Ignoring among-site heterogeneity:
    if (length(sub$Q()) > 1) {
        Q <- Reduce(`+`, sub$Q()) / length(sub$Q())
    } else Q <- sub$Q()[[1]]

    # Fill and check the position column in `x$mats()`
    mats <- fill_coal_mat_pos(x$mats(), chrom_sizes)

    haplotypes_ptr <- add_ssites_cpp(reference$ptr(),
                                   mats,
                                   Q,
                                   sub$pi_tcag(),
                                   ins$rates(),
                                   del$rates(),
                                   n_threads,
                                   show_progress)

    return(haplotypes_ptr)

}


#' Create haplotypes from VCF file
#'
#'
#' @noRd
#'
to_hap_set__haps_vcf_info <- function(x, reference, sub, ins, del, epsilon,
                                     n_threads, show_progress) {

    haplotypes_ptr <- read_vcf_cpp(reference$ptr(), x$fn(), x$print_names())

    return(haplotypes_ptr)

}


#' Create haplotypes from phylogenetic tree(s).
#'
#'
#' @noRd
#'
to_hap_set__haps_phylo_info <- function(x, reference, sub, ins, del, epsilon,
                                       n_threads, show_progress) {

    phy <- x$phylo()

    n_chroms <- as.integer(reference$n_chroms())

    if (length(phy) == 1 && n_chroms != 1) phy <- rep(phy, n_chroms)

    if (length(phy) !=  n_chroms) {
        stop("\nIn function `haps_phylo`, you must provide information for 1 tree ",
             "or a tree for each reference genome chromosome. ",
             "It appears you need to re-run `haps_phylo` before attempting to ",
             "run `create_haplotypes` again.")
    }

    trees_info <- phylo_to_info_list(phy, reference)

    hap_set_ptr <- trees_to_hap_set(trees_info, reference, sub, ins, del, epsilon,
                                    n_threads, show_progress)

    return(hap_set_ptr)

}


#' Create haplotypes from theta parameter.
#'
#'
#' @noRd
#'
to_hap_set__haps_theta_info <- function(x,
                                       reference,
                                       sub, ins, del, epsilon,
                                       n_threads, show_progress) {

    phy <- x$phylo()
    theta <- x$theta()

    n_haps <- length(phy$tip.label)
    n_chroms <- reference$n_chroms()

    # From here:
    # https://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-quantitative-
    # genomics-fall-2005/study-materials/hstnotes.pdf

    # Calculating L from theta:
    # E(L) = 4 * N * a; a = sum(1 / (1:(n_tips-1)))
    a <- sum(1 / (1:(n_haps-1)))
    # theta = 4 * N * mu
    # ------------*
    # Calculating mu:
    # ------------*
    # Indel rates (same for each nucleotide):
    indel <- sum(ins$rates() * 0.25) + sum(del$rates() * 0.25)
    # Average substution rate for each nucleotide (goes across Gammas):
    avg_subs <- -1 * colMeans(do.call(rbind, lapply(sub$Q(), diag)))
    # Average mutation rate among all nucleotides:
    mu <- sum({avg_subs + indel} * sub$pi_tcag())
    # ------------*
    # So if we know theta and mu (and since theta = 4 * N * mu), then...
    # ------------*
    L <- theta / mu * a
    # Now rescale to have total branch length of `L`:
    phy$edge.length <- phy$edge.length / sum(phy$edge.length) * L


    phy <- rep(list(phy), n_chroms)

    trees_info <- phylo_to_info_list(phy, reference)

    hap_set_ptr <- trees_to_hap_set(trees_info, reference, sub, ins, del, epsilon,
                                    n_threads, show_progress)

    return(hap_set_ptr)

}


#' Create haplotypes from gene trees.
#'
#'
#' @noRd
#'
to_hap_set__haps_gtrees_info <- function(x, reference, sub, ins, del, epsilon,
                                        n_threads, show_progress) {

    trees_info <- gtrees_to_info_list(x$trees(), reference)

    hap_set_ptr <- trees_to_hap_set(trees_info, reference, sub, ins, del, epsilon,
                                    n_threads, show_progress)

    return(hap_set_ptr)

}



# ====================================================================================`
# ====================================================================================`

# * CREATE_VARIANTS * -----

# ====================================================================================`
# ====================================================================================`

# doc start ----
#' Create haplotypes from a reference genome.
#'
#' Uses one of multiple methods to create variant haplotypes from a reference genome.
#' See \code{\link{haps_functions}} for the methods available.
#'
#'
#'
#'
#' @param reference A \code{ref_genome} object from which to generate haplotypes.
#'     This argument is required.
#' @param haps_info Output from one of the \code{\link{haps_functions}}.
#'     These functions organize higher-level information for use here.
#'     See \code{\link{haps_functions}} for brief descriptions and links to each method.
#'     If this argument is `NULL`, all arguments other than `reference` are ignored,
#'     and an empty `haplotypes` object with no haplotypes is returned.
#'     This is designed for use when you'd like to add mutations manually.
#'     If you create a blank `haplotypes` object, you can use its `add_haps` method
#'     to add haplotypes manually.
#' @param sub Output from one of the \code{\link{sub_models}} functions that organizes
#'     information for the substitution models.
#'     See \code{\link{sub_models}} for more information on these models and
#'     their required parameters.
#'     This argument is ignored if you are using a VCF file to create haplotypes.
#'     Passing `NULL` to this argument results in no substitutions.
#'     Defaults to `NULL`.
#' @param ins Output from the \code{\link{indels}} function that specifies rates
#'     of insertions by length.
#'     This argument is ignored if you are using a VCF file to create haplotypes.
#'     Passing `NULL` to this argument results in no insertions.
#'     Defaults to `NULL`.
#' @param del Output from the \code{\link{indels}} function that specifies rates
#'     of deletions by length.
#'     This argument is ignored if you are using a VCF file to create haplotypes.
#'     Passing `NULL` to this argument results in no deletions.
#'     Defaults to `NULL`.
#' @param epsilon Error control parameter for the "tau-leaping" approximation to
#'     the Doob–Gillespie algorithm, as used for the indel portion of the simulations.
#'     Smaller values result in a closer approximation.
#'     Larger values are less exact but faster.
#'     Values must be `>= 0` and `< 1`.
#'     For more information on the approximation, see Cao et al. (2006) and
#'     Wieder et al. (2011), listed below.
#'     If `epsilon` is `0`, then it reverts to the exact Doob–Gillespie algorithm.
#'     Defaults to `0.03`.
#' @param n_threads Number of threads to use for parallel processing.
#'     This argument is ignored if OpenMP is not enabled.
#'     Threads are spread across chromosomes, so it
#'     doesn't make sense to supply more threads than chromosomes in the reference genome.
#'     Defaults to `1`.
#' @param show_progress Boolean for whether to show a progress bar during processing.
#'     Defaults to `FALSE`.
#'
#'
#' @export
#'
#' @return A \code{\link{haplotypes}} object.
#'
#'
#' @references
#' Cao, Y., D. T. Gillespie, and L. R. Petzold. 2006. Efficient step size
#' selection for the tau-leaping simulation method.
#' \emph{The Journal of Chemical Physics} \strong{124}(4): 044109.
#'
#' Doob, J. L. 1942. Topics in the theory of markoff chains.
#' \emph{Transactions of the American Mathematical Society} \strong{52}(1): 37–64.
#'
#' Gillespie, D. T. 1976. A general method for numerically simulating the stochastic time
#' evolution of coupled chemical reactions. \emph{Journal of Computational Physics}
#' \strong{22}(4): 403–434.
#'
#' Wieder, N., R. H. Fink, and F. von Wegner. 2011. Exact and approximate stochastic
#' simulation of intracellular calcium dynamics.
#' \emph{Journal of Biomedicine and Biotechnology} \strong{2011}: 572492.
#'
#'
#'
#' @examples
#' r <- create_genome(10, 1000)
#' v_phylo <- create_haplotypes(r, haps_phylo(ape::rcoal(5)), sub_JC69(0.1))
#' v_theta <- create_haplotypes(r, haps_theta(0.001, 5), sub_K80(0.1, 0.2))
#'
# doc end ----
create_haplotypes <- function(reference,
                            haps_info,
                            sub = NULL,
                            ins = NULL,
                            del = NULL,
                            epsilon = 0.03,
                            n_threads = 1,
                            show_progress = FALSE) {

    # `haps_info` classes:
    vic <- list(phylo = c("phylo", "gtrees", "theta"),
                              non = c("ssites", "vcf"))
    vic <- lapply(vic, function(x) paste0("haps_", x, "_info"))

    # ---------*
    # --- check types ----
    # ---------*

    if (!inherits(reference, "ref_genome")) {
        err_msg("create_haplotypes", "reference", "a \"ref_genome\" object")
    }

    # Make empty `haplotypes` object, ignoring everything other than `reference` argument:
    if (is.null(haps_info)) {
        haplotypes_ptr <- make_hap_set(reference$ptr(), 0)
        hap_obj <- haplotypes$new(haplotypes_ptr, reference$ptr())
        return(hap_obj)
    }


    if (!inherits(reference$ptr(), "externalptr")) {
        err_msg("create_haplotypes", "reference", "a \"ref_genome\" object with a `ptr`",
                "method that returns an object of class \"externalptr\".",
                "Restart by reading a FASTA file or by simulating a genome.")
    }
    if (!inherits(haps_info, do.call(c, vic))) {
        err_msg("create_haplotypes", "haps_info", "NULL or one of the following classes:",
                paste(sprintf("\"%s\"", do.call(c, vic)), collapse = ", "))
    }
    # Check that sub, ins, or del info was passed if a non-VCF method is desired:
    vcf <- inherits(haps_info, vic$non[grepl("vcf", vic$non)])
    if (!vcf && is.null(sub) && is.null(ins) && is.null(del)) {
        stop("\nFor the `create_haplotypes` function in jackalope, ",
             "if you are using a haplotype-creation method other than a VCF file, ",
             "you must provide input to the `sub`, `ins`, or `del` argument.",
             call. = FALSE)
    }

    # -----------------*
    # Do checks and organize molecular-evolution info
    # (or `NULL` if `sub` was not provided):
    # -----------------*
    if (!is.null(sub) && !inherits(sub, "sub_info")) {
        err_msg("create_haplotypes", "sub", "NULL or a \"sub_info\" object")
    }
    if (!is.null(ins) && !inherits(ins, "indel_info")) {
        err_msg("create_haplotypes", "ins", "NULL or a \"indel_info\" object")
    }
    if (!is.null(del) && !inherits(del, "indel_info")) {
        err_msg("create_haplotypes", "del", "NULL or a \"indel_info\" object")
    }
    if (!single_number(epsilon, 0) || epsilon >= 1) {
        err_msg("create_haplotypes", "epsilon", "a single number >= 0 and < 1")
    }
    # If sub is NULL and it's not a vcf file method, convert sub to rate-0 matrix:
    if (is.null(sub) && !vcf) sub <- sub_JC69(0)

    # Below will turn `NULL` into indel_info object with `numeric(0)` as `rates` method:
    if (is.null(ins)) ins <- indel_info$new(numeric(0))
    if (is.null(del)) del <- indel_info$new(numeric(0))


    if (!single_integer(n_threads, .min = 1)) {
        err_msg("create_haplotypes", "n_threads", "a single integer >= 1")
    }
    if (!is_type(show_progress, "logical", 1)) {
        err_msg("create_haplotypes", "show_progress", "a single logical")
    }

    # `to_hap_set` is a method defined for each class of input for `haps_info`
    haplotypes_ptr <- to_hap_set(x = haps_info,
                               reference = reference,
                               sub = sub,
                               ins = ins,
                               del = del,
                               epsilon = epsilon,
                               n_threads = n_threads,
                               show_progress = show_progress)

    hap_obj <- haplotypes$new(haplotypes_ptr, reference$ptr())

    return(hap_obj)

}




