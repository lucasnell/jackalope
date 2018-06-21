context("Testing phylogenetic evolution accuracy")

library(gemino)

set.seed(1953993132)

# Create pointer to C++ sampler object
sampler <- gemino:::make_sampler(sub_params = list(pi_tcag = rep(0.25, 4),
                                                   alpha_1 = 1, alpha_2 = 1.5, beta = 2),
                                 indel_params = list(xi = 1, psi = 1,
                                                     rel_insertion_rates = exp(-1:-10),
                                                     rel_deletion_rates = exp(-1:-10)),
                                 model = "TN93", chunk_size = 100)

# Random sequences
seqs <- rando_seqs(10, 1e3)

# Phylogenetic tree:
tree <- ape::rcoal(5)
tree$edge.length <- tree$edge.length * 0.1

ordered_tip_labels <- sort(tree$tip.label)

# Using this gamma matrix sets all gamma values to 1
gamma_mat <- cbind(1000, 1)

phylo_sim_fun <- function() {
    # Pointer to VarSet object
    vars <- gemino:::make_vars(seqs, 5)
    n_muts <- as.list(0:9)
    for (i in 0:9) {
        n_muts[[i+1]] <-
            gemino:::test_phylo(
                vars,
                sampler,
                seq_ind = i,
                branch_lens = tree$edge.length,
                edges = tree$edge,
                tip_labels = tree$tip.label,
                ordered_tip_labels = ordered_tip_labels,
                gamma_mat = gamma_mat
            )
    }
    return(do.call(rbind, n_muts))
}

# Expected proportions of mutations at each edge:
expected <- tree$edge.length / sum(tree$edge.length)

# Realized proportions:
phylo_sims <- lapply(1:10, function(i) phylo_sim_fun())
phylo_sims <- do.call(rbind, phylo_sims)
phylo_sims <- t(apply(phylo_sims, 1, function(x) x / sum(x)))

test_that(paste("mutation counts on phylogeny edges are not significantly",
                "different from expectations"), {
    for (i in 1:ncol(phylo_sims)) {
        x <- phylo_sims[,i]
        p <- t.test(x, mu = expected[i])$p.value
        expect_gt(p, 0.05)
    }
})



