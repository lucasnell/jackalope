context("Testing phylogenetic evolution accuracy")


library(testthat)
library(gemino)

set.seed(1953993132)

# Create pointer to C++ sampler object
sampler <- gemino:::make_sampler(sub_params = list(pi_tcag = rep(0.25, 4),
                                                   alpha_1 = 1, alpha_2 = 1.5, beta = 2),
                                 indel_params = list(xi = 1, psi = 1,
                                                     rel_insertion_rates = exp(-1:-10),
                                                     rel_deletion_rates = exp(-1:-10)),
                                 model = "TN93", chunk_size = 100)

# Rates for each nucleotide (q vector in Yang (2006))
rates <- c(1 * 0.25 + 2 * (0.25+0.25) + 0.25 * 1,
           1 * 0.25 + 2 * (0.25+0.25) + 0.25 * 1,
           1.5 * 0.25 + 2 * (0.25+0.25) + 0.25 * 1,
           1.5 * 0.25 + 2 * (0.25+0.25) + 0.25 * 1)

# Random sequences
seqs <- rando_seqs(100, 1e3)

# Phylogenetic tree:
tree <- ape::rcoal(5)
tree$edge.length <- tree$edge.length * 0.1

ordered_tip_labels <- sort(tree$tip.label)

# Using this gamma matrix sets all gamma values to 1
gamma_mat <- cbind(1000, 1)


# Pointer to VarSet object
vars <- gemino:::make_vars(seqs, 5)

# Expected proportions of mutations at each edge:
expected <- tree$edge.length / sum(tree$edge.length)

# Realized proportions:
phylo_sims <- as.list(0:(length(seqs) - 1))
set.seed(546085105)
for (i in 0:(length(seqs) - 1)) {
    if (i < 50) {
        phylo_sims[[i+1]] <-
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
    } else {
        phylo_sims[[i+1]] <-
            gemino:::test_phylo(
                vars,
                sampler,
                seq_ind = i,
                branch_lens = tree$edge.length,
                edges = tree$edge,
                tip_labels = tree$tip.label,
                ordered_tip_labels = ordered_tip_labels,
                gamma_mat = gamma_mat,
                recombination = TRUE,
                start = 100,
                end = 499
            )
    }
}
phylo_sims <- do.call(rbind, phylo_sims)
phylo_sims <- t(apply(phylo_sims, 1, function(x) x / sum(x)))


test_that(paste("mutation counts on phylogeny edges are not significantly",
                "different from expectations"), {
    pvals <- sapply(1:ncol(phylo_sims), function(i) {
        x <- phylo_sims[,i]
        p <- t.test(x, mu = expected[i])$p.value
        return(p)
    })
    p <- pchisq(-2 * sum(log(pvals)), df = 2 * length(pvals))
    expect_gt(p, 0.05)
})



# List of data frames describing the mutations for each sequence mutated within a range:
mut_list <- lapply(0:4, function(var_ind) {
    df_ <- gemino:::see_mutations(vars, var_ind)
    return(df_[df_$seq >= 50,])
})
mut_df <- do.call(rbind, mut_list)

test_that("Ranged mutations do not occur outside specified range.", {
    expect_identical(range(mut_df$old_pos), c(100, 499))
})


proper_dels <- logical(5 * 50)
i <- 1
for (var_ind in 0:4) {
    for (seq_ind in 50:99) {
        mut_df_ <- mut_df[mut_df$var == var_ind & mut_df$seq == seq_ind,]
        max_ <- max(mut_df_$new_pos)
        max_possible <- 499 + sum(mut_df_$size_mod)

        if (max_ > max_possible) {
            if (max_ - max_possible == 1 & tail(mut_df_$size_mod, 1) < 0) {
                del_size <- abs(tail(mut_df_$size_mod, 1))
                del_start <- tail(mut_df_$old_pos, 1)
                if ((del_start + del_size - 1) == 499) {
                    proper_dels[i] <- TRUE
                } else proper_dels[i] <- FALSE
            } else proper_dels[i] <- FALSE
        } else proper_dels[i] <- TRUE
        i <- i + 1
    }
}


test_that("Ranged mutations deal with deletions at the end of sequences properly.", {
    expect_true(all(proper_dels))
})








# ------------------------------------
# ------------------------------------
# Testing that rate method works
# ------------------------------------
# ------------------------------------

# Full sequences:
var_seqs <- lapply(0:(length(tree$tip.label)-1),
                   function(v_) gemino:::see_vg(vs_ = vars, v = v_))

# R function to get expected rate to compare against C++ version:
get_seq_rate <- function(seq, rates, start, end) {
    bases <- c("T", "C", "A", "G")
    vec <- strsplit(substr(seq, start, end), "")[[1]]
    rates_ <- sapply(vec, function(x) rates[which(bases == x)])
    return(sum(rates_))
}

# Compare R version to C++ version by using random start and end points
compare_rates <- function(var_ind, seq_ind) {
    nchars <- nchar(var_seqs[[var_ind]][seq_ind])
    start <- sample.int(as.integer(nchars / 2), 1)
    end <- sample.int(nchars - start, 1) + start
    rate_cpp <- gemino:::test_rate(start = start-1, end = end-1,
                                   var_ind = var_ind-1, seq_ind = seq_ind-1,
                                   var_set_sexp = vars, sampler_sexp = sampler)
    rate_r <- get_seq_rate(var_seqs[[var_ind]][seq_ind], rates, start, end)
    if (rate_cpp != rate_r) cat(sprintf("%i, %i, %i, %i\n", var_ind, seq_ind, start, end))
    return(cbind(rate_cpp, rate_r))
}


mat <- matrix(0, 500, 2)
i <- 1
for (var_ind in 1:5) {
    for (seq_ind in 1:100) {
        mat[i,] <- compare_rates(var_ind, seq_ind)
        i <- i + 1
    }
}

test_that("Rate method produces accurate results", {
    expect_identical(mat[,1], mat[,2])
})

