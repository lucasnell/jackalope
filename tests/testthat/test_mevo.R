context("Testing molecular evolution accuracy")

library(gemino)

set.seed(1087437799)
# Construct sequence with known number of each nucleotide
seq <- c(rep("T", 0.25e6), rep("C", 0.25e6), rep("A", 0.25e6), rep("G", 0.25e6))
seq <- sample(seq)
seq <- paste(seq, collapse = "")

# Make pointer to RefGenome object based on `seq`
ref <- gemino:::make_ref_genome(seq)

# Set molecular evolution parameters
N_ <- 1e3
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

gm <- matrix(1L, 1e3, 2)
gm[,1] <- seq(1e3, 1e6, 1e3)
set.seed(1087437801)
gm[,2] <- rgamma(1e3, shape = 1, rate = 1)


# Make a matrix have ten columns
make_ten <- function(x) cbind(x, matrix(0, nrow(x), 10 - ncol(x)))


test_samp <- function() {
    vars <- gemino:::make_var_set(ref, 1)
    gemino:::test_sampling(var_set_ = vars,
                           N = N_,
                           pi_tcag = pi_tcag_,
                           alpha_1 = alpha_1_,
                           alpha_2 = alpha_2_,
                           beta = beta_,
                           xi = xi_,
                           psi = psi_,
                           rel_insertion_rates = exp(-1:-10),
                           rel_deletion_rates = exp(-1:-10),
                           gamma_mat = gm,
                           chunk_size = 100,
                           display_progress = FALSE)
    mut_exam <- gemino:::examine_mutations(vars, 0, 0)
    mut_exam$ins <- make_ten(mut_exam$ins)
    mut_exam$del <- make_ten(mut_exam$del)
    return(mut_exam)
}




nreps <- 100

set.seed(1255574685)
muts <- lapply(1:nreps, function(i) test_samp())



# Checking the indel accuracy:

all_ins <- do.call(rbind, lapply(1:nreps, function(i) setNames(colSums(muts[[i]]$ins),
                                                               1:10)))
all_del <- do.call(rbind, lapply(1:nreps, function(i) setNames(colSums(muts[[i]]$del),
                                                               1:10)))
# Function to calculate the expected number of insertions or deletions of a given size:
indel_p <- function(ss) (N_ * 0.5 * xi_ / sum(q)) * (exp(-ss) / sum(exp(-1:-10)))



test_that("molecular evolution produces correct proportion of indel sizes", {
    expect_true(t.test(rowSums(all_ins) / rowSums(all_del), mu = psi_)$p.value > 0.05)
    ins_pvals <- sapply(1:10, function(i) {
        # If it sums to > 0, then we'll do simple t.test
        if (sum(all_ins[,i]) > 0) {
            p <- t.test(all_ins[,i], mu = indel_p(i))$p.value
        # If none were sampled, do binomial test for probability
        } else {
            p <- binom.test(0, N_ * nreps, p = indel_p(i) / N_)$p.value
        }
        return(p)
    })
    # Combined P value for all items in matrix (using Fisher's method):
    p <- pchisq(-2 * sum(log(ins_pvals)), df = 2 * 10)
    expect_true(p > 0.05, info = sprintf("insertion size: %i", i))

    # Same for deletions
    del_pvals <- sapply(1:10, function(i) {
        # If it sums to > 0, then we'll do simple t.test
        if (sum(all_del[,i]) > 0) {
            p <- t.test(all_del[,i], mu = indel_p(i))$p.value
            # If none were sampled, do binomial test for probability
        } else {
            p <- binom.test(0, N_ * nreps, p = indel_p(i) / N_)$p.value
        }
        return(p)
    })
    # Combined P value for all items in matrix (using Fisher's method):
    p <- pchisq(-2 * sum(log(del_pvals)), df = 2 * 10)
    expect_true(p > 0.05, info = sprintf("deletion size: %i", i))
})




# Checking the substition accuracy:


sub_df <- lapply(1:nreps,
                 function(i) {
                     data.frame(r = i,
                                count = as.integer(muts[[i]]$sub),
                                from = rep(1:4, 4),
                                to = rep(1:4, each = 4))
                 })
sub_df <- do.call(rbind, sub_df)
sub_df <- sub_df[sub_df$from != sub_df$to,]

expected_sub <- function(i, j) {
    N_ * ((q[i] - 0.25 * xi_) / sum(q)) * (Q[i, j] / sum(Q[i,]))
}

# Observed mean counts per run
Q_obs <- matrix(0, 4, 4)
for (i in 1:4) {
    for (j in 1:4) {
        if (i == j) next
        Q_obs[i, j] <- mean(sub_df$count[sub_df$from == i & sub_df$to == j])
    }
}
# Expected
Q_exp <- t(sapply(1:4, function(i) sapply(1:4, function(j) expected_sub(i, j))))

pvals <- numeric(12)
k <- 0
for (i in 1:4) {
    for (j in 1:4) {
        if (i == j) next
        k <- k + 1
        pvals[k] <- t.test(sub_df$count[sub_df$from == i & sub_df$to == j],
                           mu = Q_exp[i,j])$p.value
    }
}

# Combined P value for all items in matrix (using Fisher's method):
p <- pchisq(-2 * sum(log(pvals)), df = 2 * k)


test_that("molecular evolution produces correct substitution matrix", {
    # I'm not just using alpha in case there's positive correlation
    threshold <- 0.05 * (k+1) / (2*k)
    expect_true(p > threshold)
})




# Now looking at positions

pos_df <- lapply(1:nreps,
                 function(i) {
                     data.frame(rep = i,
                                gamma = gm[,2],
                                pos = gm[,1],
                                count = gemino:::table_gammas(gm[,1], muts[[i]]$pos))
                 })
pos_df <- do.call(rbind, pos_df)

pos_df <- lapply(1:nrow(gm), function(i) {
    pos_ <- gm[i,1]
    gamma_ <- gm[i,2]
    count_ <- mean(pos_df$count[pos_df$pos == pos_])
    data.frame(pos = pos_, gamma = gamma_, count = count_)
})
pos_df <- do.call(rbind, pos_df)


# This regression coefficient should approximately correspond to one
test_that("molecular evolution selects mutation regions according to Gamma values", {
    gamma_coef <- coef(lm(count ~ gamma, data = pos_df))[['gamma']]
    expect_true(gamma_coef > 0.9 & gamma_coef < 1.1)
})





