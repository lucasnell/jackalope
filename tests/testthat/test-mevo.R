
context("Testing molecular evolution accuracy")

# library(gemino)
# library(testthat)



# simulate ----

set.seed(1087437799)

# Construct sequences with known number of each nucleotide
n_seqs <- 20
seqs <- rep(list(c(rep("T", 0.25e3), rep("C", 0.25e3),
                   rep("A", 0.25e3), rep("G", 0.25e3))),
            n_seqs)
seqs <- lapply(seqs, sample)
seqs <- sapply(seqs, paste, collapse = "")

# Make ref_genome object from a pointer to a RefGenome object based on `seqs`
ref <- ref_genome$new(gemino:::make_ref_genome(seqs))

# Set all needed molecular evolution parameters inside an environment (using TN93 method)
pars <- new.env()
with(pars, {
    pi_tcag = c(0.1, 0.2, 0.3, 0.4)
    alpha_1 = 0.25
    alpha_2 = 0.35
    beta = 0.5
    # For indels:
    rates = c(0.2, 0.3)
    M = c(8L, 10L)
    a = c(0.1, 0.5)
    # For site variability:
    shape = 0.5
    region_size = 100
    seq_len = nchar(seqs[1])
    ends = seq(region_size, seq_len, region_size)
    mats = replicate(n_seqs,
                     cbind(ends, rgamma(length(ends), shape = shape, rate = shape)),
                     simplify = FALSE)
})


mevo_ <- with(pars, {
    make_mevo(ref,
              sub = list(model = "TN93", alpha_1 = alpha_1, alpha_2 = alpha_2,
                         beta = beta, pi_tcag = pi_tcag),
              ins = list(rate = rates[1], max_length = M[1], a = a[1]),
              del = list(rate = rates[2], max_length = M[2], a = a[2]),
              site_var = list(mats = pars$mats),
              chunk_size = 0)
})

n_vars <- 10

# Phylogenetic trees:
trees <- lapply(1:n_seqs,
                function(i) {
                    tree <- ape::rcoal(n_vars)
                    # Don't want this to take too long:
                    tree$edge.length <- tree$edge.length * 0.1
                    return(tree)
                })

var_set <- create_variants(ref, method = "phylo", method_info = trees,
                           mevo_obj = mevo_)





# =================================================================*
# =================================================================*

# examine output ----

# =================================================================*
# =================================================================*


# Function to make a matrix have N columns
make_N <- function(x, N) cbind(x, matrix(0, nrow(x), N - ncol(x)))

mutation_exam <- function(v, s) {
    mut_exam <- gemino:::examine_mutations(var_set, var_ind = v, seq_ind = s)
    mut_exam$ins <- setNames(colSums(make_N(mut_exam$ins, pars$M[1])), 1:pars$M[1])
    mut_exam$del <- setNames(colSums(make_N(mut_exam$del, pars$M[2])), 1:pars$M[2])
    return(mut_exam)
}


exams <- lapply(0:(n_seqs-1),
                function(s) lapply(0:(n_vars-1),
                                    gemino:::examine_mutations,
                                    var_set_ = var_set$genomes,
                                    seq_ind = s))


# --------------------*
# indels ----
# --------------------*

# Function to calculate the expected proportion of mutations that are
# insertions or deletions of a given size:
indel_p <- function(u, j) {
    stopifnot(!missing(j))
    # < indel of type j > / < all mutations >
    ijp <- pars$rates[j] / sum(mevo_$q)
    # < indel of type j and size u > / < indel of type j >
    M_ <- pars$M[j]
    a_ <- pars$a[j]
    u_ <- 1:M_
    all_ps <- {(u_ * M_) / (M_ - u_ + 1)}^(-a_)
    piju <- all_ps[u] / sum(all_ps)
    return(ijp * piju)
}

# Sample for only one variant to avoid the issue of closely related
# individuals having many of the same mutations:
rand_vars <- sample.int(n_vars, n_seqs, replace = TRUE)

# Coercing total # mutations (# trials) to a single vector
total_muts <- sapply(1:n_seqs,
                       function(s) {
                           v <- rand_vars[s]
                           length(exams[[s]][[v]]$pos)
                       })

# Columns are different sizes:
insertions <- sapply(1:pars$M[1], function(u) {
    sapply(1:n_seqs,
             function(s) {
                 v <- rand_vars[s]
                 if (ncol(exams[[s]][[v]]$ins) > pars$M[1]) {
                     M <- exams[[s]][[v]]$ins[,1:pars$M[1]]
                 } else M <- exams[[s]][[v]]$ins
                 sum(make_N(M, pars$M[1])[,u])
             })
})
deletions <- sapply(1:pars$M[2], function(u) {
    sapply(1:n_seqs,
             function(s) {
                 v <- rand_vars[s]
                 if (ncol(exams[[s]][[v]]$del) > pars$M[2]) {
                     M <- exams[[s]][[v]]$del[,1:pars$M[2]]
                 } else M <- exams[[s]][[v]]$del
                 sum(make_N(M, pars$M[2])[,u])
             })
})



test_that("molecular evolution produces correct proportion of indel sizes", {
    # Combining p-values using Fisher's method:
    fisher_method <- function(pvals) {
        chisq <- -2 * sum(log(pvals))
        df_ <- 2 * length(pvals)
        p <- pchisq(chisq, df = df_, lower.tail = FALSE)
        return(p)
    }
    ins_pvals <- sapply(1:pars$M[1], function(u) {
        binom.test(sum(insertions[,u]), sum(total_muts), p = indel_p(u, 1))$p.value
    })
    # Combined P value for all items in matrix (using Fisher's method):
    p <- fisher_method(ins_pvals)
    alpha_hat <- 0.05 / length(ins_pvals)  # <-- Bonferroni correction
    expect_gt(p, alpha_hat, label = "insertion sizes")

    # Same for deletions
    del_pvals <- sapply(1:pars$M[2], function(u) {
        binom.test(sum(deletions[,u]), sum(total_muts), p = indel_p(u, 2))$p.value
    })
    # Combined P value for all items in matrix (using Fisher's method):
    p <- fisher_method(del_pvals)
    alpha_hat <- 0.05 / length(del_pvals)  # <-- Bonferroni correction
    expect_gt(p, alpha_hat, label = "deletion sizes")
})



# --------------------*
# substitutions ----
# --------------------*

sub_df <- do.call(
    rbind,
    lapply(1:n_seqs,
           function(s) {
               # Sample for only one variant to avoid the issue of closely related
               # individuals having many of the same mutations:
               v <- rand_vars[s]
               data.frame(seq = s, var = v,
                          from = rep(1:4, 4),
                          to = rep(1:4, each = 4),
                          n = sum(exams[[s]][[v]]$sub),
                          count = as.integer(exams[[s]][[v]]$sub))
           }))
sub_df <- sub_df[sub_df$from != sub_df$to,]

# Expected proportion
expected_sub <- function(i, j) {
    Q <- mevo_$Q
    q <- rowSums(Q)
    (q[i] / sum(q)) * (Q[i, j] / sum(Q[i,]))
}

# Observed counts per variant--sequence combo
Q_obs <- matrix(0, 4, 4)
# Expected:
Q_exp <- matrix(0, 4, 4)
for (i in 1:4) {
    for (j in 1:4) {
        if (i == j) next
        df_ <- sub_df[sub_df$from == i & sub_df$to == j,]
        Q_obs[i, j] <- sum(df_$count)
        Q_exp[i, j] <- expected_sub(i, j) * sum(df_$n)
    }
}

test_that("molecular evolution produces correct substitution matrix", {
    x <- c(Q_obs[lower.tri(Q_obs)], Q_obs[upper.tri(Q_obs)])
    y <- c(Q_exp[lower.tri(Q_exp)], Q_exp[upper.tri(Q_exp)])
    p <- t.test(x, y, paired = TRUE)$p.value
    expect_gt(p, 0.05)
})




# --------------------*
# gammas ----
# --------------------*

pos_df <- do.call(
    rbind,
    lapply(1:n_seqs,
           function(s) {
               # Sample for only one variant to avoid the issue of closely related
               # individuals having many of the same mutations:
               v <- rand_vars[s]
               exam_ <- exams[[s]][[v]]
               gamma_mat <- mevo_$gamma_mats[[s]]
               # Expected # mutations per gamma region:
               exp_mut = length(exam_$pos) / nrow(gamma_mat)
               df_ <- data.frame(seq = s, var = v,
                                 gamma = gamma_mat[,2],
                                 pos = gamma_mat[,1],
                                 count = gemino:::table_gammas(gamma_mat[,1], exam_$pos))
               df_$count <- df_$count / exp_mut
               return(df_)
           }))


# This regression coefficient should approximately correspond to one
test_that("molecular evolution selects mutation regions according to Gamma values", {
    mod <- lm(count ~ gamma + 0, data = pos_df)
    gamma_coef <- coef(mod)[['gamma']]
    expect_true(gamma_coef > 0.5 & gamma_coef < 1.5)
})






# --------------------*
# phylogeny ----
# --------------------*


mutations <- lapply(0:(n_vars-1), gemino:::view_mutations, var_set_ = var_set$genomes)


# Function to return descendent tips from a given node:
get_descendant_tips <- function(tree, node) {
    # Function to get descendent nodes and tips from a given node
    # From http://blog.phytools.org/2012/01/function-to-get-descendant-node-numbers.html
    # Slight tweeks for personal style
    getDescendants <- function(tree, node, curr = NULL){
        if(is.null(curr)) curr <- vector()
        daughters <- tree$edge[which(tree$edge[,1] == node), 2]
        curr <- c(curr, daughters)
        w <- which(daughters >= length(tree$tip))
        if (length(w) > 0) {
            for (i in 1:length(w)) {
                curr <- getDescendants(tree, daughters[w[i]], curr)
            }
        }
        return(curr)
    }
    tips_ <- getDescendants(tree, node)
    tips_ <- tips_[tips_ <= length(tree$tip.label)]
    return(tips_)
}

compare_mutations <- function(df1, df2) {

    df1 <- df1[, c("seq", "old_pos", "size_mod", "nucleos")]
    df2 <- df2[, c("seq", "old_pos", "size_mod", "nucleos")]

    df1_ <- apply(df1, 1, function(xx) gsub(" ", "", paste(xx, collapse = "_")))
    df2_ <- apply(df2, 1, function(xx) gsub(" ", "", paste(xx, collapse = "_")))

    n_shared_muts <- sum(df2_ %in% df1_)

    return(as.integer(n_shared_muts))

}


shared_muts <- rep(list(matrix(0L, n_vars, n_vars)), n_seqs)
for (s in 1:n_seqs) {
    for (i in 1:(n_vars-1)) {
        for (j in (i+1):n_vars) {
            df1 <- mutations[[i]]
            df2 <- mutations[[j]]
            sm <- compare_mutations(df1[df1$seq == (s-1), ], df2[df2$seq == (s-1), ])
            shared_muts[[s]][i,j] <- sm
            shared_muts[[s]][j,i] <- sm
        }
    }
    diag(shared_muts[[s]]) <- sapply(mutations, function(x) sum(x$seq==(s-1)))
}; rm(s, i, j, df1, df2, sm)

# "Mutational" distances between all pairs:
mut_dists <- rep(list(matrix(0L, n_vars, n_vars)), n_seqs)
for (s in 1:n_seqs) {
    for (i in 1:(n_vars-1)) {
        for (j in (i+1):n_vars) {
            # Total non-shared mutations for each:
            mdi <- shared_muts[[s]][i,i] - shared_muts[[s]][i,j]
            mdj <- shared_muts[[s]][j,j] - shared_muts[[s]][i,j]
            mut_dists[[s]][i,j] <- mdi
            mut_dists[[s]][j,i] <- mdj
        }
    }
}; rm(s, i, j, mdi, mdj)


# Phylogenetic distances between all pairs:
phylo_dists <- lapply(trees, function(tree) {
    pd <- ape::cophenetic.phylo(tree)
    rownames(pd) <- colnames(pd) <- NULL
    # Dividing by two to get branch lengths between MRCA:
    pd <- pd / 2
    return(pd)
})



# Average mutation rate for ancestral sequences times total bp in each sequence:
mu_ <- mean(mevo_$q) * nchar(seqs[1])

# Expected mutational distances:
expected <- lapply(phylo_dists, function(pd) mu_ * pd)


s = 18

mapply(expected, mut_dists,
       FUN = function(x, y) mean((x - y) / x, na.rm = TRUE))

mean((expected[[s]] - mut_dists[[s]]) / expected[[s]], na.rm = TRUE)



mean(sapply(mut_dists, function(x) mean((expected - x) / expected,
                                        na.rm = TRUE)))


test_that("Simulations conform to phylogenetic tree", {
    expect_lt(mean(sapply(mut_dists, function(x) mean((expected - x) / expected,
                                                      na.rm = TRUE))), expected = 0.1)
})





# =================================================================*
# =================================================================*

# rate method ----

# =================================================================*
# =================================================================*

# Full sequences:
var_seqs <- lapply(0:(n_vars-1),
                   function(v) gemino:::view_var_genome(var_set_ = var_set, var_ind = v))

# R function to get expected rate to compare against C++ version:
get_seq_rate <- function(seq, rates, start, end, gamma_mat_) {
    bases <- c("T", "C", "A", "G")
    vec <- strsplit(substr(seq, start, end), "")[[1]]
    rates_ <- sapply(vec, function(x) rates[which(bases == x)])
    gammas_ <- sapply(start:end, function(i) gamma_mat_[,2][which(gamma_mat_[,1] >= i)[1]])
    return(sum(rates_ * gammas_))
}

# Compare R version to C++ version by using random start and end points
compare_rates <- function(var_ind, seq_ind, rates, incr = 1e2) {
    nchars <- nchar(var_seqs[[var_ind]][seq_ind])
    if (nchars <= (2 * incr)) stop("set incr to a lower value")
    gamma_ends <- rev(seq(nchars, incr, -incr))
    gamma_mat_new <- cbind(gamma_ends,
                           rgamma(length(gamma_ends), shape = 1, rate = 1))
    start <- sample.int(as.integer(nchars / 2), 1)
    end <- sample.int(nchars - start, 1) + start
    rate_cpp <- gemino:::test_rate(start = start-1, end = end-1,
                                   var_ind = var_ind-1, seq_ind = seq_ind-1,
                                   var_set_ = var_set, sampler_ = sampler,
                                   gamma_mat_ = gamma_mat_new)
    rate_r <- get_seq_rate(var_seqs[[var_ind]][seq_ind], rates, start, end, gamma_mat_new)
    if (!all.equal(rate_cpp, rate_r)) cat(sprintf("%i, %i, %i, %i\n", var_ind, seq_ind, start, end))
    return(cbind(rate_cpp, rate_r))
}


mat <- matrix(0, n_vars * n_seqs, 2)
i <- 1
for (var_ind in 1:n_vars) {
    for (seq_ind in 1:n_seqs) {
        mat[i,] <- compare_rates(var_ind, seq_ind, q)
        i <- i + 1
    }
}

test_that("Rate method produces accurate results", {
    expect_equal(mat[,1], mat[,2])
})




