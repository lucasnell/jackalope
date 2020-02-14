
# library(jackalope)
# library(testthat)


context("Testing making of mevo (molecular evolution info) object")


set.seed(4616515)


# Set all needed molecular evolution parameters inside an environment
pars <- new.env()
with(pars, {
    # For reference genome:
    n_chroms <- 10
    chrom_len <- 1000
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
})

# Create reference genome
ref <- with(pars, create_genome(n_chroms, chrom_len))




# ==============================*
# __sub mats__ ----
# ==============================*

# *  JC69 ----
M <- sub_JC69(pars$lambda)


test_that("proper output for JC69 model", {
    Q <- matrix(pars$lambda, 4, 4)
    diag(Q) <- 0
    diag(Q) <- -1 * rowSums(Q)
    Q <- Q / abs(sum(diag(Q) * rep(0.25, 4)))  # Scale to overall mutation rate of 1
    expect_equal(M$Q()[[1]], Q)
    expect_equal(M$pi_tcag(), rep(0.25, 4))
})

# Same but no scaling:
M <- sub_JC69(pars$lambda, mu = NULL)
test_that("proper output for JC69 model", {
    Q <- matrix(pars$lambda, 4, 4)
    diag(Q) <- 0
    diag(Q) <- -1 * rowSums(Q)
    expect_equal(M$Q()[[1]], Q)
    expect_equal(M$pi_tcag(), rep(0.25, 4))
})

# *  K80 ----
M <- with(pars, {sub_K80(alpha = alpha, beta = beta)})

test_that("proper output for K80 model", {
    Q <- matrix(pars$beta, 4, 4)
    Q[2,1] <- Q[1,2] <- Q[4,3] <- Q[3,4] <- pars$alpha
    diag(Q) <- 0
    diag(Q) <- -1 * rowSums(Q)
    Q <- Q / abs(sum(diag(Q) * rep(0.25, 4)))  # Scale to overall mutation rate of 1
    expect_equal(M$Q()[[1]], Q)
    expect_equal(M$pi_tcag(), rep(0.25, 4))
})


# *  F81 ----
M <- with(pars, {sub_F81(pi_tcag = pi_tcag)})

test_that("proper output for F81 model", {
    Q <- matrix(rep(pars$pi_tcag, each = 4), 4, 4)
    diag(Q) <- 0
    diag(Q) <- -1 * rowSums(Q)
    Q <- Q / abs(sum(diag(Q) * pars$pi_tcag))  # Scale to overall mutation rate of 1
    expect_equal(M$Q()[[1]], Q)
    expect_equal(M$pi_tcag(), pars$pi_tcag)
})


# *  HKY85 ----
M <- with(pars, {sub_HKY85(alpha = alpha, beta = beta,
                                           pi_tcag = pi_tcag)})

test_that("proper output for HKY85 model", {
    Q <- matrix(pars$beta, 4, 4)
    Q[2,1] <- Q[1,2] <- Q[4,3] <- Q[3,4] <- pars$alpha
    for (i in 1:4) Q[,i] <- Q[,i] * pars$pi_tcag[i]
    diag(Q) <- 0
    diag(Q) <- -1 * rowSums(Q)
    Q <- Q / abs(sum(diag(Q) * pars$pi_tcag))  # Scale to overall mutation rate of 1
    expect_equal(M$Q()[[1]], Q)
    expect_equal(M$pi_tcag(), pars$pi_tcag)
})

# *  TN93 ----
M <- with(pars, {sub_TN93(alpha_1 = alpha_1,
                                           alpha_2 = alpha_2, beta = beta,
                                           pi_tcag = pi_tcag)})

test_that("proper output for TN93 model", {
    Q <- matrix(pars$beta, 4, 4)
    Q[2,1] <- Q[1,2] <- pars$alpha_1
    Q[4,3] <- Q[3,4] <- pars$alpha_2
    for (i in 1:4) Q[,i] <- Q[,i] * pars$pi_tcag[i]
    diag(Q) <- 0
    diag(Q) <- -1 * rowSums(Q)
    Q <- Q / abs(sum(diag(Q) * pars$pi_tcag))  # Scale to overall mutation rate of 1
    expect_equal(M$Q()[[1]], Q)
    expect_equal(M$pi_tcag(), pars$pi_tcag)
})

# *  F84 ----
M <- with(pars, {sub_F84(beta = beta, kappa = kappa,
                                           pi_tcag = pi_tcag)})

test_that("proper output for F84 model", {
    alpha_1 <- (1 + pars$kappa / sum(pars$pi_tcag[1:2])) * pars$beta
    alpha_2 <- (1 + pars$kappa / sum(pars$pi_tcag[3:4])) * pars$beta
    Q <- matrix(pars$beta, 4, 4)
    Q[2,1] <- Q[1,2] <- alpha_1
    Q[4,3] <- Q[3,4] <- alpha_2
    for (i in 1:4) Q[,i] <- Q[,i] * pars$pi_tcag[i]
    diag(Q) <- 0
    diag(Q) <- -1 * rowSums(Q)
    Q <- Q / abs(sum(diag(Q) * pars$pi_tcag))  # Scale to overall mutation rate of 1
    expect_equal(M$Q()[[1]], Q)
    expect_equal(M$pi_tcag(), pars$pi_tcag)
})

# *  GTR ----
M <- with(pars, {sub_GTR(pi_tcag = pi_tcag, abcdef = abcdef)})

test_that("proper output for GTR model", {
    Q <- matrix(0, 4, 4)
    Q[lower.tri(Q)] <- pars$abcdef
    Q <- Q + t(Q)
    for (i in 1:4) Q[,i] <- Q[,i] * pars$pi_tcag[i]
    diag(Q) <- -1 * rowSums(Q)
    Q <- Q / abs(sum(diag(Q) * pars$pi_tcag))  # Scale to overall mutation rate of 1
    expect_equal(M$Q()[[1]], Q)
    expect_equal(M$pi_tcag(), pars$pi_tcag)

    # Special case when all are zeros:
    Msp <- sub_GTR(pi_tcag = pars$pi_tcag, abcdef = pars$abcdef * 0)
    expect_equal(Msp$Q()[[1]], matrix(0, 4, 4))
})

# *  UNREST ----
M <- with(pars, {sub_UNREST(Q = Q)})

test_that("proper output for UNREST model", {
    Q <- pars$Q
    diag(Q) <- 0
    diag(Q) <- -1 * rowSums(Q)
    eig <- eigen(t(Q))
    eig_vals <- abs(eig$values)
    eig_vecs <- Re(eig$vectors)
    i <- which(eig_vals == min(eig_vals))
    left_vec <- eig_vecs[,i]
    sumlv <- sum(left_vec)
    pi_tcag <- left_vec / sumlv
    Q <- Q / abs(sum(diag(Q) * pi_tcag))  # Scale to overall mutation rate of 1
    expect_equal(M$Q()[[1]], Q)
    expect_equal(M$pi_tcag(), pi_tcag)
})



# * errors ----

test_that("proper errors for sub models", {

    # pi_tcag = NULL, alpha_1 = NULL, alpha_2 = NULL,
    # beta = NULL, gamma_shape = NULL, gamma_k = NULL,
    # invariant = NULL, lambda = NULL, alpha = NULL,
    # kappa = NULL, abcdef = NULL, Q = NULL

    expect_error(sub_JC69(lambda = "lambda"), "`lambda` must be a single number >= 0")
    expect_error(sub_K80(alpha = function(x) x, beta = pars$beta),
                 "`alpha` must be a single number >= 0")
    expect_error(sub_TN93(alpha_1 = "alpha_1", alpha_2 = pars$alpha_2, beta = pars$beta,
                          pi_tcag = pars$pi_tcag),
                 "`alpha_1` must be a single number >= 0")
    expect_error(sub_TN93(alpha_1 = pars$alpha_1, alpha_2 = -1, beta = pars$beta,
                          pi_tcag = pars$pi_tcag),
                 "`alpha_2` must be a single number >= 0")
    expect_error(sub_TN93(alpha_1 = pars$alpha_1, alpha_2 = pars$alpha_2, beta = "beta",
                          pi_tcag = pars$pi_tcag),
                 "`beta` must be a single number >= 0")
    expect_error(sub_F84(beta = pars$beta, kappa = -1, pi_tcag = pars$pi_tcag),
                 "`kappa` must be a single number >= 0")

    err <- paste("`pi_tcag` must be a length-4 numeric vector where at least one",
                 "number is > 0 and all are >= 0")
    expect_error(sub_TN93(alpha_1 = pars$alpha_1, alpha_2 = pars$alpha_2,
                          beta = pars$beta, pi_tcag = pars$pi_tcag[1]), err)
    expect_error(sub_TN93(alpha_1 = pars$alpha_1, alpha_2 = pars$alpha_2,
                          beta = pars$beta, pi_tcag = pars$pi_tcag[1:3]), err)
    expect_error(sub_TN93(alpha_1 = pars$alpha_1, alpha_2 = pars$alpha_2,
                          beta = pars$beta, pi_tcag = pars$pi_tcag * -1), err)
    expect_error(sub_TN93(alpha_1 = pars$alpha_1, alpha_2 = pars$alpha_2,
                          beta = pars$beta, pi_tcag = pars$pi_tcag * 0), err)


    err <- paste("`abcdef` must be a length-6 numeric vector where all numbers are >= 0")
    expect_error(sub_GTR(pi_tcag = pars$pi_tcag, abcdef = pars$abcdef * -1), err)
    expect_error(sub_GTR(pi_tcag = pars$pi_tcag, abcdef = pars$abcdef[1:3]), err)

    err <- paste("`Q` must be a 4x4 numeric matrix where all non-diagonal elements",
                 "are >= 0")
    expect_error(sub_UNREST(Q = "pars$Q * -1"), err)
    expect_error(sub_UNREST(Q = pars$Q * -1), err)
    expect_error(sub_UNREST(Q = pars$Q[1:3,]), err)
    expect_error(sub_UNREST(Q = pars$Q[2:4]), err)

})



# ==============================*
# __indel rates__ ----
# ==============================*

# *  exp(-L) ----
test_that("proper indel rates with `rate` and `max_length` inputs", {
    R <- with(pars, indels(rate = rates[1], max_length = M[1]))
    ins <- exp(-1 * 1:pars$M[1])
    ins <- (ins / sum(ins)) * pars$rates[1]
    expect_equal(R$rates(), ins)

    R <- with(pars, indels(rate = rates[2], max_length = M[2]))
    del <- exp(-1 * 1:pars$M[2])
    del <- (del / sum(del)) * pars$rates[2]
    expect_equal(R$rates(), del)
})

# *  Lavalette ----
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

    R <- with(pars, indels(rate = rates[1], max_length = M[1], a = a[1]))
    expect_equal(R$rates(), ins)
    R <- with(pars, indels(rate = rates[2], max_length = M[2], a = a[2]))
    expect_equal(R$rates(), del)

})

# * custom ----
test_that("proper indel rates with `rate` and `rel_rates` inputs", {
    ins <- pars$rel_rates[[1]]
    ins <- (ins / sum(ins)) * pars$rates[1]
    del <- pars$rel_rates[[2]]
    del <- (del / sum(del)) * pars$rates[2]
    R <- with(pars, indels(rate = rates[1], rel_rates = rel_rates[[1]]))
    expect_equal(R$rates(), ins)
    R <- with(pars, indels(rate = rates[2], rel_rates = rel_rates[[2]]))
    expect_equal(R$rates(), del)
})




# ==============================*
# __site var.__ ----
# ==============================*


M <- with(pars, {sub_TN93(alpha_1 = alpha_1,
                          alpha_2 = alpha_2, beta = beta,
                          pi_tcag = pi_tcag,
                          gamma_shape = 1)})

# *  rate matrices ----
test_that("proper substitution rate matrices for discrete Gamma model", {

    Q <- with(pars, sub_TN93(alpha_1 = alpha_1, alpha_2 = alpha_2,
                             beta = beta, pi_tcag = pi_tcag))$Q()[[1]]

    # Each in M$Q() should be simply result of multiplying Q by gamma category
    gammas <- lapply(M$Q(), function(x) as.numeric(x / Q))
    # Make sure they're all roughly the same within each matrix:
    expect_lt(max(sapply(gammas, function(x) max(diff(abs(x))))), 1e-10)

    # Mean for each cell in matrices within M$Q() list should be equal to its cell in Q
    expect_equal(Reduce(`+`, M$Q()) / length(M$Q()), Q)

    # Default # categories:
    expect_length(M$Q(), 5)

    expect_equal(M$invariant(), 0.0)

})



# *  invariants ----
test_that("invariant sites are produced appropriately", {
    # Because `ref` is evenly sized we should get exact output:
    M <- sub_JC69(0.1, gamma_shape = 1, invariant = 0.5)
    expect_equal(M$invariant(), 0.5)
})




# *  proper errors ----
test_that("throws proper errors when inputting nonsensible gamma input", {

    err <- "`gamma_shape` must be NULL or a single number > 0"
    expect_error(sub_JC69(0.1, gamma_shape = 0), regexp = err)
    expect_error(sub_JC69(0.1, gamma_shape = -1), regexp = err)
    expect_error(sub_JC69(0.1, gamma_shape = c(1, 2)), regexp = err)
    expect_error(sub_JC69(0.1, gamma_shape = "0"), regexp = err)

    err <- "`gamma_k` must be a single integer in range \\[2, 255\\]"
    expect_error(sub_JC69(0.1, gamma_shape = 1, gamma_k = "-1"), regexp = err)
    expect_error(sub_JC69(0.1, gamma_shape = 1, gamma_k = -1), regexp = err)
    expect_error(sub_JC69(0.1, gamma_shape = 1, gamma_k = 1), regexp = err)
    expect_error(sub_JC69(0.1, gamma_shape = 1, gamma_k = 256), regexp = err)


    err <- "`invariant` must be a single number >= 0 and < 1"
    expect_error(sub_JC69(0.1, gamma_shape = 1, invariant = -1), regexp = err)
    expect_error(sub_JC69(0.1, gamma_shape = 1, invariant = c(0.1, 0.5)), regexp = err)
    expect_error(sub_JC69(0.1, gamma_shape = 1, invariant = 1), regexp = err)
})

