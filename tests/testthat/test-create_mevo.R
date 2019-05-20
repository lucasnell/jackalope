
context("Testing making of mevo (molecular evolution info) object")

# library(jackalope)
# library(testthat)



set.seed(4616515)


# Set all needed molecular evolution parameters inside an environment
pars <- new.env()
with(pars, {
    # For reference genome:
    n_seqs <- 10
    seq_len <- 1000
    # Molecular evolution:
    lambda = 0.1
    alpha = 0.25
    beta = 0.5
    pi_tcag = c(0.1, 0.2, 0.3, 0.4)
    alpha_1 = 0.25
    alpha_2 = 0.35
    kappa = 0.75
    abcdef = 1:6 * 0.1
    Q = matrix(1:16, 4, 4)
    # For indels:
    rates = c(0.2, 0.3)
    M = c(8L, 10L)
    a = c(0.1, 0.5)
    rel_rates = list(1:10 * 0.1, 8:1 * 0.2)
    # For site variability:
    shape = 0.5
    region_size = 10
    ends = seq(region_size, seq_len, region_size)
    mats = replicate(n_seqs,
                     cbind(ends, rgamma(length(ends), shape = shape, rate = shape)),
                     simplify = FALSE)
    # Gamma matrices that should throw errors:
    mats_err1 = replicate(n_seqs,
                          cbind(ends, rgamma(length(ends), shape = shape, rate = shape)),
                          simplify = TRUE)
    mats_err2 = mats
    mats_err2[[2]] = mats_err2[[2]][-nrow(mats_err2[[2]]),]
    mats_err3 = mats
    mats_err3[[3]][,1] = mats_err3[[3]][,1] + 0.1
    mats_err4 = mats
    mats_err4[[4]][1,1] = mats_err4[[4]][2,1]
})

# This will be used a lot
create_mevo <- function(reference, sub,
                        ins = NULL,
                        del = NULL,
                        gamma_mats = NULL) {
    jackalope:::create_mevo(reference, sub, ins, del, gamma_mats, 10)
}

# Create reference genome
ref <- with(pars, create_genome(n_seqs, seq_len))




# ==============================*
# __sub mats__ ----
# ==============================*

# *  JC69 ----
M <- create_mevo(ref, sub = sub_JC69(pars$lambda))

# These only need to be checked once:
test_that("proper output common to all models with no site variability or indels", {

    X <- lapply(ref$sizes(),
                function(x) {
                    z <- cbind(x, 1)
                    colnames(z) <- NULL
                    return(z)
                })

    expect_equal(M$gamma_mats, X, check.attributes = FALSE)
    expect_equal(M$insertion_rates, numeric(0))
    expect_equal(M$deletion_rates, numeric(0))
})

test_that("proper output for JC69 model", {
    Q <- matrix(pars$lambda, 4, 4)
    diag(Q) <- 0
    expect_equal(M$Q, Q)
    expect_equal(M$pi_tcag, rep(0.25, 4))
    expect_equal(M$q(), rowSums(Q))
})

# *  K80 ----
M <- with(pars, {create_mevo(ref, sub = sub_K80(alpha = alpha, beta = beta))})

test_that("proper output for K80 model", {
    Q <- matrix(pars$beta, 4, 4)
    Q[2,1] <- Q[1,2] <- Q[4,3] <- Q[3,4] <- pars$alpha
    diag(Q) <- 0
    expect_equal(M$Q, Q)
    expect_equal(M$pi_tcag, rep(0.25, 4))
    expect_equal(M$q(), rowSums(Q))
})


# *  F81 ----
M <- with(pars, {create_mevo(ref, sub_F81(pi_tcag = pi_tcag))})

test_that("proper output for F81 model", {
    Q <- matrix(rep(pars$pi_tcag, each = 4), 4, 4)
    diag(Q) <- 0
    expect_equal(M$Q, Q)
    expect_equal(M$pi_tcag, pars$pi_tcag)
    expect_equal(M$q(), rowSums(Q))
})


# *  HKY85 ----
M <- with(pars, {create_mevo(ref, sub_HKY85(alpha = alpha, beta = beta,
                                           pi_tcag = pi_tcag))})

test_that("proper output for HKY85 model", {
    Q <- matrix(pars$beta, 4, 4)
    Q[2,1] <- Q[1,2] <- Q[4,3] <- Q[3,4] <- pars$alpha
    for (i in 1:4) Q[,i] <- Q[,i] * pars$pi_tcag[i]
    diag(Q) <- 0
    expect_equal(M$Q, Q)
    expect_equal(M$pi_tcag, pars$pi_tcag)
    expect_equal(M$q(), rowSums(Q))
})

# *  TN93 ----
M <- with(pars, {create_mevo(ref, sub_TN93(alpha_1 = alpha_1,
                                           alpha_2 = alpha_2, beta = beta,
                                           pi_tcag = pi_tcag))})

test_that("proper output for TN93 model", {
    Q <- matrix(pars$beta, 4, 4)
    Q[2,1] <- Q[1,2] <- pars$alpha_1
    Q[4,3] <- Q[3,4] <- pars$alpha_2
    for (i in 1:4) Q[,i] <- Q[,i] * pars$pi_tcag[i]
    diag(Q) <- 0
    expect_equal(M$Q, Q)
    expect_equal(M$pi_tcag, pars$pi_tcag)
    expect_equal(M$q(), rowSums(Q))
})

# *  F84 ----
M <- with(pars, {create_mevo(ref, sub_F84(beta = beta, kappa = kappa,
                                           pi_tcag = pi_tcag))})

test_that("proper output for F84 model", {
    alpha_1 <- (1 + pars$kappa / sum(pars$pi_tcag[1:2])) * pars$beta
    alpha_2 <- (1 + pars$kappa / sum(pars$pi_tcag[3:4])) * pars$beta
    Q <- matrix(pars$beta, 4, 4)
    Q[2,1] <- Q[1,2] <- alpha_1
    Q[4,3] <- Q[3,4] <- alpha_2
    for (i in 1:4) Q[,i] <- Q[,i] * pars$pi_tcag[i]
    diag(Q) <- 0
    expect_equal(M$Q, Q)
    expect_equal(M$pi_tcag, pars$pi_tcag)
    expect_equal(M$q(), rowSums(Q))
})

# *  GTR ----
M <- with(pars, {create_mevo(ref, sub_GTR(pi_tcag = pi_tcag,
                                           abcdef = abcdef))})

test_that("proper output for GTR model", {
    Q <- matrix(0, 4, 4)
    Q[lower.tri(Q)] <- pars$abcdef
    Q <- Q + t(Q)
    for (i in 1:4) Q[,i] <- Q[,i] * pars$pi_tcag[i]
    expect_equal(M$Q, Q)
    expect_equal(M$pi_tcag, pars$pi_tcag)
    expect_equal(M$q(), rowSums(Q))
})

# *  UNREST ----
M <- with(pars, {create_mevo(ref, sub_UNREST(Q = Q))})

test_that("proper output for UNREST model", {
    Q <- pars$Q
    diag(Q) <- 0  # <-- important for eigen step
    diag(Q) <- -1 * rowSums(Q)  # <-- important for eigen step
    eig <- eigen(t(Q))
    eig_vals <- abs(eig$values)
    eig_vecs <- Re(eig$vectors)
    i <- which(eig_vals == min(eig_vals))
    left_vec <- eig_vecs[,i]
    sumlv <- sum(left_vec)
    pi_tcag <- left_vec / sumlv
    diag(Q) <- 0  # <-- change back to compare to M$Q
    expect_equal(M$Q, Q)
    expect_equal(M$pi_tcag, pi_tcag)
    expect_equal(M$q(), rowSums(Q))
})


# ==============================*
# __indel rates__ ----
# ==============================*

# *  exp(-L) ----
M <- with(pars, create_mevo(ref, sub_JC69(lambda = lambda),
                          ins = indels(rate = rates[1], max_length = M[1]),
                          del = indels(rate = rates[2], max_length = M[2])))
test_that("proper indel rates with `rate` and `max_length` inputs", {
    ins <- exp(-1 * 1:pars$M[1])
    ins <- (ins / sum(ins)) * pars$rates[1]
    del <- exp(-1 * 1:pars$M[2])
    del <- (del / sum(del)) * pars$rates[2]
    expect_equal(M$insertion_rates, ins)
    expect_equal(M$deletion_rates, del)
})

# *  Lavalette ----
M <- with(pars, create_mevo(ref, sub_JC69(lambda = lambda),
                          ins = indels(rate = rates[1], max_length = M[1], a = a[1]),
                          del = indels(rate = rates[2], max_length = M[2], a = a[2])))
test_that("proper indel rates with `rate`, `max_length`, and `a` inputs", {
    u <- 1:pars$M[1]
    M_ <- pars$M[1]
    a <- pars$a[1]
    ins <- {(u * M_) / (M_ - u + 1)}^(-a)
    ins <- (ins / sum(ins)) * pars$rates[1]
    u <- 1:pars$M[2]
    M_ <- pars$M[2]
    a <- pars$a[2]
    del <- {(u * M_) / (M_ - u + 1)}^(-a)
    del <- (del / sum(del)) * pars$rates[2]
    expect_equal(M$insertion_rates, ins)
    expect_equal(M$deletion_rates, del)
})

# * custom ----
M <- with(pars, create_mevo(ref, sub_JC69(lambda = lambda),
                          ins = indels(rate = rates[1], rel_rates = rel_rates[[1]]),
                          del = indels(rate = rates[2], rel_rates = rel_rates[[2]])))
test_that("proper indel rates with `rate` and `rel_rates` inputs", {
    ins <- pars$rel_rates[[1]]
    ins <- (ins / sum(ins)) * pars$rates[1]
    del <- pars$rel_rates[[2]]
    del <- (del / sum(del)) * pars$rates[2]
    expect_equal(M$insertion_rates, ins)
    expect_equal(M$deletion_rates, del)
})



test_that("proper output for JC69 model when indels are included", {
    indel <- sum(pars$rates * 0.25)
    Q <- matrix(pars$lambda, 4, 4)
    diag(Q) <- 0
    expect_equal(M$Q, Q)
    expect_equal(M$pi_tcag, rep(0.25, 4))
    expect_equal(M$q(), rowSums(Q) + indel)
})




# ==============================*
# __site var.__ ----
# ==============================*

dir <- tempdir(check = TRUE)

M <- site_var(ref, shape = pars$shape, region_size = pars$region_size,
              out_prefix = paste0(dir, "/mevo"))

# *  generate ----
test_that("proper gamma distance values with `shape` and `region_size` inputs", {
    # Testing end points:
    expect_equal(sapply(M, function(x) x[,1]),
                 matrix(rep(seq(pars$region_size, pars$seq_len, pars$region_size),
                            pars$n_seqs),
                        ceiling(pars$seq_len / pars$region_size), pars$n_seqs))
    # Testing mean:
    expect_equal(mean(sapply(M, function(x) mean(x[,2]))), 1)
    # Testing SD:
    G <- do.call(c, lapply(M, function(x) x[,2]))
    s <- sd(G)                 # observed SD
    s0 <- sqrt(1 / pars$shape)      # expected SD
    rel_diff <- abs((s - s0) / s0)
    expect_lte(rel_diff, 0.25)  # was never >0.25 in 1000 sims
})



test_that("BED file of gamma values is correct", {

    bed_df <- utils::read.table(paste0(dir, "/mevo.bed"), sep = "\t")

    gmat_df <- do.call(rbind, lapply(1:length(M),
                                     function(i) {
                                         cbind(nm = ref$names()[i],
                                               as.data.frame(M[[i]]))
                                         }))

    expect_equal(bed_df[,1], gmat_df[,1])
    expect_equal(bed_df[,3], gmat_df[,2])
    expect_equal(bed_df[,5], gmat_df[,3])

})


# *  custom ----

M <- site_var(ref, mats = pars$mats)

test_that("proper gamma distance values with `mats` inputs", {
    MM <- M
    class(MM) <- "list"
    expect_identical(MM, pars$mats)
})


# *  proper errors ----
test_that("throws proper errors when inputting an incorrect `site_var$mats` input", {
    expect_error(site_var(ref, mats = pars$mats_err1),
                 regexp = "argument `mats` must be NULL or a list of matrices")
    expect_error(site_var(ref, mats = pars$mats_err2),
                 regexp = paste("all matrices need to have a maximum end point",
                                "\\(in the first column\\) equal to the size of",
                                "the associated sequence"))
    expect_error(site_var(ref, mats = pars$mats_err3),
                 regexp = "all matrices should contain only whole numbers as end points")
    expect_error(site_var(ref, mats = pars$mats_err4),
                 regexp = "all matrices should contain no duplicate end points")
})



# *  invariants ----
test_that("invariant sites are produced appropriately", {
    # Because `ref` is evenly sized we should get exact output:
    M <- site_var(ref, shape = pars$shape, region_size = pars$region_size,
                  invariant = 0.50)
    expect_identical(sapply(M, function(x) mean(x[,2] == 0)),
                     rep(0.5, length(M)))
    expect_equal(mean(sapply(M, function(x) mean(x[,2]))), 1.0)
})
