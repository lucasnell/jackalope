context("Testing molecular evolution accuracy")

library(gemino)


# Combining p-values using Fisher's method:
fisher_method <- function(pvals) {
    chisq <- -2 * sum(log(pvals))
    df_ <- 2 * length(pvals)
    p <- pchisq(chisq, df = df_, lower.tail = FALSE)
    return(p)
}



# =================================================================
# =================================================================

#  Do the simulations including a phylogeny and evolutionary model

# =================================================================
# =================================================================



set.seed(1087437799)


# Construct sequences with known number of each nucleotide
n_seqs <- 20
seqs <- rep(list(c(rep("T", 0.25e3), rep("C", 0.25e3),
                   rep("A", 0.25e3), rep("G", 0.25e3))),
            n_seqs)
seqs <- lapply(seqs, sample)
seqs <- sapply(seqs, paste, collapse = "")

# Make pointer to RefGenome object based on `seqs`
ref <- gemino:::make_ref_genome(seqs)

# Set molecular evolution parameters
pi_tcag_ = 0.1*1:4
alpha_1_ = 1
alpha_2_ = 2
beta_ = 0.5
xi_ = 0.5
psi_ = 1

Q <- matrix(beta_, 4, 4)
Q[2,1] <- Q[1,2] <- alpha_1_
Q[4,3] <- Q[3,4] <- alpha_2_
for (i in 1:4) Q[,i] <- Q[,i] * pi_tcag_[i]
diag(Q) <- 0

q <- rowSums(Q) + 0.25 * xi_

gamma_mat <- matrix(1L, 10L, 2)
gamma_mat[,1] <- seq(1e2, 1e3, 1e2)
gamma_mat[,2] <- rgamma(nrow(gamma_mat), shape = 1, rate = 1)
# Repeating the above gamma matrix for half of all sequences, bc we're
# processing that many at a time:
gamma_mats <- rep(list(gamma_mat), n_seqs / 2)

# Create pointer to C++ sampler object
sampler <- gemino:::make_sampler(sub_params = list(pi_tcag = pi_tcag_,
                                                   alpha_1 = alpha_1_,
                                                   alpha_2 = alpha_2_,
                                                   beta = beta_),
                                 indel_params = list(xi = xi_, psi = psi_,
                                                     rel_insertion_rates = exp(-1:-10),
                                                     rel_deletion_rates = exp(-1:-10)),
                                 model = "TN93", chunk_size = 100)


n_vars <- 10

# Phylogenetic tree:
tree <- ape::rcoal(n_vars)
tree$edge.length <- tree$edge.length * 0.1
# Expected proportions of mutations at each edge:
expected <- tree$edge.length / sum(tree$edge.length)
ordered_tip_labels <- sort(tree$tip.label)



# Function to make a matrix have ten columns
make_ten <- function(x) cbind(x, matrix(0, nrow(x), 10 - ncol(x)))


var_set <- gemino:::make_var_set(ref, n_vars)

edge_muts <- as.list(0:(n_seqs - 1))
edge_muts[1:(n_seqs/2)] <- gemino:::test_mevo(
    var_set_ = var_set,
    sampler_base_ = sampler,
    seq_inds = 0:(n_seqs/2 - 1),
    branch_lens = tree$edge.length,
    edges = tree$edge,
    tip_labels = tree$tip.label,
    ordered_tip_labels = ordered_tip_labels,
    gamma_mats = gamma_mats
)
# Do the second half within a range (as would be done w recombination)
edge_muts[(n_seqs/2+1):n_seqs] <- gemino:::test_mevo(
    var_set_ = var_set,
    sampler_base_ = sampler,
    seq_inds = (n_seqs/2):(n_seqs-1),
    branch_lens = tree$edge.length,
    edges = tree$edge,
    tip_labels = tree$tip.label,
    ordered_tip_labels = ordered_tip_labels,
    gamma_mats = gamma_mats,
    recombination = TRUE,
    start = 100,
    end = 499
)





# =================================================================
# =================================================================

# Checking evolutionary model output

# =================================================================
# =================================================================



mutation_exam <- function(v, s) {
    mut_exam <- gemino:::examine_mutations(var_set, var_ind = v, seq_ind = s)
    mut_exam$ins <- setNames(colSums(make_ten(mut_exam$ins)), 1:10)
    mut_exam$del <- setNames(colSums(make_ten(mut_exam$del)), 1:10)
    return(mut_exam)
}

exams <- lapply(0:(n_seqs-1),
                function(s) lapply(0:(n_vars-1),
                                    gemino:::examine_mutations,
                                    var_set_ = var_set,
                                    seq_ind = s))


# --------------------
# Checking the indel accuracy:
# --------------------

# Function to calculate the expected proportion of mutations that are
# insertions or deletions of a given size:
indel_p <- function(u) {
    (0.5 * xi_ / sum(q)) * (exp(-u) / sum(exp(-1:-10)))
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
insertions <- sapply(1:10, function(u) {
    sapply(1:n_seqs,
             function(s) {
                 v <- rand_vars[s]
                 sum(make_ten(exams[[s]][[v]]$ins)[,u])
             })
})
deletions <- sapply(1:10, function(u) {
    sapply(1:n_seqs,
             function(s) {
                 v <- rand_vars[s]
                 sum(make_ten(exams[[s]][[v]]$del)[,u])
             })
})



test_that("molecular evolution produces correct proportion of indel sizes", {
    expect_true(t.test(rowSums(insertions), rowSums(deletions))$p.value > 0.05)
    ins_pvals <- sapply(1:10, function(u) {
        binom.test(sum(insertions[,u]), sum(total_muts), p = indel_p(u))$p.value
    })
    # Combined P value for all items in matrix (using Fisher's method):
    p <- fisher_method(ins_pvals)
    expect_true(p > 0.05, info = "insertion sizes")

    # Same for deletions
    del_pvals <- sapply(1:10, function(u) {
        binom.test(sum(deletions[,u]), sum(total_muts), p = indel_p(u))$p.value
    })
    # Combined P value for all items in matrix (using Fisher's method):
    p <- fisher_method(del_pvals)
    expect_true(p > 0.05, info = "deletion sizes")
})



# --------------------
# Checking the substition accuracy:
# --------------------


sub_df <- do.call(
    rbind,
    lapply(1:n_seqs,
           function(s) {
               # Sample for only one variant to avoid the issue of closely related
               # individuals having many of the same mutations:
               v <- rand_vars[s]
               data.frame(seq = s, var = v,
                          n = length(exams[[s]][[v]]$pos),
                          count = as.integer(exams[[s]][[v]]$sub),
                          from = rep(1:4, 4),
                          to = rep(1:4, each = 4))
           }))
sub_df <- sub_df[sub_df$from != sub_df$to,]

# Expected proportion
expected_sub <- function(i, j) {
    ((q[i] - 0.25 * xi_) / sum(q)) * (Q[i, j] / sum(Q[i,]))
}

# Observed counts per variant--sequence combo
Q_obs <- matrix(0, 4, 4)
# Number of mutations per variant--sequence combo
n_obs <- matrix(0, 4, 4)
# Expected:
Q_exp <- matrix(0, 4, 4)
for (i in 1:4) {
    for (j in 1:4) {
        if (i == j) next
        df_ <- sub_df[sub_df$from == i & sub_df$to == j,]
        Q_obs[i, j] <- sum(df_$count)
        n_obs[i, j] <- sum(df_$n)
        Q_exp[i, j] <- expected_sub(i, j)
    }
}
# Observed proportions:
# apply(Q_obs, c(1,2), mean) / apply(n_obs, c(1,2), mean)
# Expected:
# Q_exp

sub_pvals <- numeric(12)
k <- 1
for (i in 1:4) {
    for (j in 1:4) {
        if (i == j) next
        sub_pvals[k] <- binom.test(Q_obs[i, j], n_obs[i, j], p = Q_exp[i,j])$p.value
        k <- k + 1
    }
}; rm(k, i, j)

# Combined P value for all items in matrix (using Fisher's method):
p <- fisher_method(sub_pvals)


test_that("molecular evolution produces correct substitution matrix", {
    # I'm not just using alpha in case there's positive correlation
    k <- length(sub_pvals)
    threshold <- 0.05 * (k+1) / (2*k)
    expect_true(p > threshold)
})




# --------------------
# Gamma region accuracy
# --------------------

pos_df <- do.call(
    rbind,
    lapply(1:n_seqs,
           function(s) {
               # Sample for only one variant to avoid the issue of closely related
               # individuals having many of the same mutations:
               v <- rand_vars[s]
               exam_ <- exams[[s]][[v]]
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
    mod <- lm(count ~ gamma, data = pos_df)
    gamma_coef <- coef(mod)[['gamma']]
    expect_true(gamma_coef > 0.75 & gamma_coef < 1.25)
})







# =================================================================
# =================================================================

# Testing phylogenetic aspects of simulations

# =================================================================
# =================================================================

# Converting # mutations per edge to matrix
edge_muts <- do.call(rbind, edge_muts)

test_that(paste("mutation counts on phylogeny edges are not significantly",
                "different from expectations"), {
                    pvals <- sapply(1:ncol(edge_muts), function(i) {
                        x <- edge_muts[,i]
                        p <- binom.test(sum(x), sum(edge_muts), p = expected[i])$p.value
                        return(p)
                    })
                    p <- fisher_method(pvals)
                    expect_gt(p, 0.05)
                })



# List of data frames describing the mutations for each sequence mutated within a range:
mut_list <- lapply(0:(n_vars-1), function(var_ind) {
    df_ <- gemino:::see_mutations(var_set, var_ind)
    return(df_[df_$seq >= (n_seqs/2),])
})
mut_df <- do.call(rbind, mut_list)

test_that("Ranged mutations do not occur outside specified range.", {
    expect_gte(min(mut_df$old_pos), 100)
    expect_lte(max(mut_df$old_pos), 499)
})


proper_dels <- logical(n_vars * (n_seqs / 2))
i <- 1
for (var_ind in 0:(n_vars-1)) {
    for (seq_ind in (n_seqs/2):(n_seqs-1)) {
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






# =================================================================
# =================================================================

# Testing that rate method works

# =================================================================
# =================================================================

# Full sequences:
var_seqs <- lapply(0:(n_vars-1),
                   function(v) gemino:::see_var_genome(var_set_ = var_set, var_ind = v))

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




