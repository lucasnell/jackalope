
context("Testing molecular evolution accuracy")

# library(gemino)
# library(testthat)


boot <- function(x, B = 2000, alpha = 0.01) {
    boots <- numeric(B)
    N <- length(x)
    for (i in 1:B) {
        boots[i] <- mean(sample(x, size = N, replace = TRUE))
    }
    return(quantile(boots, c(alpha / 2, 1 - alpha / 2), names = FALSE))
}



# simulate ----

# Set all needed molecular evolution parameters inside an environment (using TN93 method)
pars <- new.env()
with(pars, {
    n_seqs <- 20
    # Molecular evolution parameters:
    pi_tcag = c(0.1, 0.2, 0.3, 0.4)
    alpha_1 = 0.25
    alpha_2 = 0.35
    beta = 0.5
    #   For indels:
    rates = c(0.2, 0.3)
    M = c(8L, 10L)
    a = c(0.1, 0.5)
    n_vars = 10
    #   For site variability:
    shape = 0.5
    region_size = 100
    seq_len = 1000
    ends = seq(region_size, seq_len, region_size)
    mats = replicate(n_seqs,
                     cbind(ends, rgamma(length(ends), shape = shape, rate = shape)),
                     simplify = FALSE)
    # Construct sequences with known number of each nucleotide
    seqs <- rep(list(c(rep("T", pi_tcag[1] * seq_len), rep("C", pi_tcag[2] * seq_len),
                       rep("A", pi_tcag[3] * seq_len), rep("G", pi_tcag[4] * seq_len))),
                n_seqs)
    seqs <- lapply(seqs, sample)
    seqs <- sapply(seqs, paste, collapse = "")
})


set.seed(1087437799)

# Make ref_genome object from a pointer to a RefGenome object based on `seqs`
ref <- with(pars, ref_genome$new(gemino:::make_ref_genome(seqs)))


mevo_ <- with(pars, {
    make_mevo(ref,
              sub = list(model = "TN93", alpha_1 = alpha_1, alpha_2 = alpha_2,
                         beta = beta, pi_tcag = pi_tcag),
              ins = list(rate = rates[1], max_length = M[1], a = a[1]),
              del = list(rate = rates[2], max_length = M[2], a = a[2]),
              site_var = list(mats = pars$mats),
              chunk_size = 100)
})


# Phylogenetic tree:
tree <- ape::rcoal(pars$n_vars)
# Don't want this to take too long:
tree$edge.length <- tree$edge.length * 0.1



var_set <- create_variants(ref, method = "phylo", method_info = tree,
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


exams <- lapply(0:(pars$n_seqs-1),
                function(s) lapply(0:(pars$n_vars-1),
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
    ijp <- pars$rates[j] / sum(mevo_$q())
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
rand_vars <- sample.int(pars$n_vars, pars$n_seqs, replace = TRUE)

# Coercing total # mutations (# trials) to a single vector
total_muts <- sapply(1:pars$n_seqs,
                       function(s) {
                           v <- rand_vars[s]
                           length(exams[[s]][[v]]$pos)
                       })

# Columns are different sizes:
insertions <- sapply(1:pars$M[1], function(u) {
    sapply(1:pars$n_seqs,
             function(s) {
                 v <- rand_vars[s]
                 if (ncol(exams[[s]][[v]]$ins) > pars$M[1]) {
                     M <- exams[[s]][[v]]$ins[,1:pars$M[1]]
                 } else M <- exams[[s]][[v]]$ins
                 sum(make_N(M, pars$M[1])[,u])
             })
})
deletions <- sapply(1:pars$M[2], function(u) {
    sapply(1:pars$n_seqs,
             function(s) {
                 v <- rand_vars[s]
                 if (ncol(exams[[s]][[v]]$del) > pars$M[2]) {
                     M <- exams[[s]][[v]]$del[,1:pars$M[2]]
                 } else M <- exams[[s]][[v]]$del
                 sum(make_N(M, pars$M[2])[,u])
             })
})



test_that("molecular evolution produces correct proportion of indel sizes", {

    ins_exp <- matrix(rep(total_muts, pars$M[1]), pars$n_seqs) *
        matrix(rep(indel_p(1:pars$M[1], 1), each = pars$n_seqs), pars$n_seqs)
    p <- chisq.test(insertions, ins_exp, simulate.p.value = TRUE, B = 1000)$p.value
    expect_gt(p, 0.01, label = "insertion sizes")

    # Same for deletions
    del_exp <- matrix(rep(total_muts, pars$M[2]), pars$n_seqs) *
        matrix(rep(indel_p(1:pars$M[2], 2), each = pars$n_seqs), pars$n_seqs)
    p <- chisq.test(deletions, del_exp, simulate.p.value = TRUE, B = 1000)$p.value
    expect_gt(p, 0.01, label = "deletion sizes")
})



# --------------------*
# substitutions ----
# --------------------*

sub_df <- do.call(
    rbind,
    lapply(1:pars$n_seqs,
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
    expect_gt(p, 0.01)
})




# --------------------*
# gammas ----
# --------------------*

pos_df <- do.call(
    rbind,
    lapply(1:pars$n_seqs,
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


mutations <- lapply(0:(pars$n_vars-1), gemino:::view_mutations, var_set_ = var_set$genomes)


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


shared_muts <- rep(list(matrix(0L, pars$n_vars, pars$n_vars)), pars$n_seqs)
for (s in 1:pars$n_seqs) {
    for (i in 1:(pars$n_vars-1)) {
        for (j in (i+1):pars$n_vars) {
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
mut_dists <- rep(list(matrix(0L, pars$n_vars, pars$n_vars)), pars$n_seqs)
for (s in 1:pars$n_seqs) {
    for (i in 1:(pars$n_vars-1)) {
        for (j in (i+1):pars$n_vars) {
            # Total non-shared mutations for each:
            mdi <- shared_muts[[s]][i,i] - shared_muts[[s]][i,j]
            mdj <- shared_muts[[s]][j,j] - shared_muts[[s]][i,j]
            mut_dists[[s]][i,j] <- mdi
            mut_dists[[s]][j,i] <- mdj
        }
    }
}; rm(s, i, j, mdi, mdj)


# Phylogenetic distances between all pairs:
phylo_dists <- ape::cophenetic.phylo(tree)
rownames(phylo_dists) <- colnames(phylo_dists) <- NULL
# Dividing by two to get branch lengths between MRCA:
phylo_dists <- phylo_dists / 2


# Average mutation rate for ancestral sequences times total bp in each sequence:
mu_ <- mevo_$mu() * pars$seq_len

# Expected mutational distances:
expected <- mu_ * phylo_dists

comp_mats <- lapply(mut_dists, function(y) (expected - y) / expected)


test_that("Simulations conform to phylogenetic tree", {
    # Bootstrap these bc they're a weird distribution
    b <- boot(sapply(comp_mats, mean, na.rm = TRUE))
    expect_lt(b[1], 0)
    expect_gt(b[2], 0)
})





# =================================================================*
# =================================================================*

# rate method ----

# =================================================================*
# =================================================================*

# Full sequences:
var_seqs <- lapply(0:(pars$n_vars-1),
                   function(v) gemino:::view_var_genome(var_set_ = var_set$genomes, var_ind = v))

# R function to get expected rate to compare against C++ version:
get_seq_rate <- function(seq, rates, start, end, gamma_mat_) {
    bases <- c("T", "C", "A", "G")
    vec <- strsplit(substr(seq, start, end), "")[[1]]
    rates_ <- sapply(vec, function(x) rates[which(bases == x)])
    gammas_ <- sapply(start:end, function(i) gamma_mat_[,2][which(gamma_mat_[,1] >= i)[1]])
    return(sum(rates_ * gammas_))
}

# Compare R version to C++ version by using random start and end points
compare_rates <- function(var_ind, seq_ind, rates, incr = 1e2, verbose = FALSE) {
    nchars <- nchar(var_seqs[[var_ind]][seq_ind])
    if (nchars <= (2 * incr)) stop("set incr to a lower value")
    gamma_ends <- rev(seq(nchars, incr, -incr))
    gamma_mat_new <- cbind(gamma_ends,
                           rgamma(length(gamma_ends), shape = 1, rate = 1))
    start <- sample.int(as.integer(nchars / 2), 1)
    end <- sample.int(nchars - start, 1) + start
    if (mevo_$chunk_size <= 0) stop("Must only be used for chunked mevo objects")
    rate_cpp <- gemino:::test_rate(start = start-1, end = end-1,
                                   var_ind = var_ind-1, seq_ind = seq_ind-1,
                                   var_set_ = var_set$genomes, sampler_ = mevo_$to_ptr(),
                                   gamma_mat_ = gamma_mat_new)
    rate_r <- get_seq_rate(var_seqs[[var_ind]][seq_ind], rates, start, end, gamma_mat_new)
    if (verbose & !isTRUE(all.equal(rate_cpp, rate_r))) {
        cat(sprintf("%i, %i, %i, %i\n", var_ind, seq_ind, start, end))
    }
    return(cbind(rate_cpp, rate_r))
}


mat <- matrix(0, pars$n_vars * pars$n_seqs, 2)
i <- 1
for (var_ind in 1:pars$n_vars) {
    for (seq_ind in 1:pars$n_seqs) {
        mat[i,] <- compare_rates(var_ind, seq_ind, rates = mevo_$q())
        i <- i + 1
    }
}

test_that("Rate method produces accurate results", {
    expect_equal(mat[,1], mat[,2])
})




