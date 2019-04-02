

library(jackalope)
library(tidyverse)
library(scrm)
library(grid)
library(ape)
source(".Rprofile")



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
    region_size = 1e3
    seq_len = 100e3
    ends = seq(region_size, seq_len, region_size)
    mats = replicate(n_seqs,
                     cbind(ends, rgamma(length(ends), shape = shape, rate = shape)),
                     simplify = FALSE)
    # Construct sequences with known number of each nucleotide
    seqs <- rep(list(c(rep("T", 0.25 * seq_len), rep("C", 0.25 * seq_len),
                       rep("A", 0.25 * seq_len), rep("G", 0.25 * seq_len))),
                n_seqs)
    seqs <- lapply(seqs, sample)
    seqs <- sapply(seqs, paste, collapse = "")
})


set.seed(1087437799)

# Make ref_genome object from a pointer to a RefGenome object based on `seqs`
ref <- with(pars, ref_genome$new(jackalope:::make_ref_genome(seqs)))



# For testing: just one individual at known sites
coal_obj <- scrm("2 20 -r 3.1 1000 -t 1000")
coal_obj$seg_sites <- map(coal_obj$seg_sites, ~ .x[1,.x[1,] > 0,drop=FALSE])




# =================================================================*
# =================================================================*

# INDELS ----

# =================================================================*
# =================================================================*


# Function to calculate the expected proportion of mutations that are
# insertions or deletions of a given size:
indel_pred <- function(u, j) {
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

# Function to make a matrix have N columns
make_N <- function(x, N) cbind(x, matrix(0, nrow(x), N - ncol(x)))

# --------------------*
# Insertions ----
# --------------------*

# Only insertions:
# ref <- create_genome(20, 100e3)
mevo_ <- make_mevo(ref,
                   sub = list(model = "TN93", alpha_1 = 0,
                              alpha_2 = 0, beta = 0,
                              pi_tcag = pars$pi_tcag),
                   ins = list(rate = pars$rates[1],
                              max_length = pars$M[1],
                              a = pars$a[1]))
var_set <- create_variants(ref, method = "coal_sites", method_info = coal_obj,
                           mevo_obj = mevo_)

insertions <- map(0:(pars$n_seqs-1),
             function(s) {
                 Z <- jackalope:::examine_mutations(var_set_ptr = var_set$genomes,
                                                 var_ind = 0, seq_ind = s) %>%
                     .[["ins"]]
                 # Anything over `pars$M[1]` is >1 mutation
                 if (ncol(Z) > pars$M[1]) Z <- Z[,1:pars$M[1],drop=FALSE]
                 # Make sure it has exactly pars$M[1] columns for binding
                 if (ncol(Z) < pars$M[1]) Z <- make_N(Z, pars$M[1])
                 return(colSums(Z))
             })%>%
    do.call(what = rbind) %>%
    colSums()


ins_df <- tibble(pred = indel_pred(1:pars$M[1], 1) * sum(insertions),
                 obs = insertions)

lims <- range(c(ins_df$obs, ins_df$pred))

ins_p <- ins_df %>%
    ggplot(aes(pred, obs)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "firebrick2") +
    scale_x_continuous("Predicted insertion count",
                       limits = lims) +
    scale_y_continuous("Observed insertion count",
                       limits = lims) +
    theme_classic()




# --------------------*
# Deletions ----
# --------------------*

# Only deletions
mevo_ <- make_mevo(ref,
                   sub = list(model = "TN93", alpha_1 = 0,
                              alpha_2 = 0, beta = 0,
                              pi_tcag = pars$pi_tcag),
                   del = list(rate = pars$rates[2],
                              max_length = pars$M[2],
                              a = pars$a[2]))
var_set <- create_variants(ref, method = "coal_sites", method_info = coal_obj,
                           mevo_obj = mevo_)

deletions <- map(0:(pars$n_seqs-1),
                  function(s) {
                      Z <- jackalope:::examine_mutations(var_set_ptr = var_set$genomes,
                                                      var_ind = 0, seq_ind = s) %>%
                          .[["del"]]
                      # Anything over `pars$M[1]` is >1 mutation
                      if (ncol(Z) > pars$M[2]) Z <- Z[,1:pars$M[2],drop=FALSE]
                      # Make sure it has exactly pars$M[1] columns for binding
                      if (ncol(Z) < pars$M[2]) Z <- make_N(Z, pars$M[2])
                      return(colSums(Z))
                  })%>%
    do.call(what = rbind) %>%
    colSums()


del_df <- tibble(pred = indel_pred(1:pars$M[2], 2) * sum(deletions),
                 obs = deletions)

lims <- range(c(del_df$obs, del_df$pred))

del_p <- del_df %>%
    ggplot(aes(pred, obs)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "firebrick2") +
    scale_x_continuous("Predicted deletion count",
                       limits = lims) +
    scale_y_continuous("Observed deletion count",
                       limits = lims) +
    theme_classic()



grid.newpage()
grid.draw(rbind(ggplotGrob(del_p), ggplotGrob(ins_p), size = "last"))




# =================================================================*
# =================================================================*

# SUBSTITUTIONS ----

# =================================================================*
# =================================================================*



# Phylogenetic tree:
tree <- ape::rcoal(2)
tree$edge.length <- tree$edge.length * 0.1

# Only substitutions
mevo_ <- make_mevo(ref,
                   sub = list(model = "TN93", alpha_1 = pars$alpha_1,
                              alpha_2 = pars$alpha_2, beta = pars$beta,
                              pi_tcag = pars$pi_tcag))
var_set <- create_variants(ref, method = "phylo", method_info = tree,
                           mevo_obj = mevo_)

substitutions <- map(0:(pars$n_seqs-1),
                     function(s) {
                         Z <- jackalope:::examine_mutations(var_set_ptr = var_set$genomes,
                                                     var_ind = 0, seq_ind = s) %>%
                         .[["sub"]]
                     return(Z)
                 }) %>%
    Reduce(f = `+`)
# Predicted values:
sub_pred <- sum(substitutions) * (mevo_$Q / sum(mevo_$Q))

sub_df <- tibble(obs = c(substitutions[lower.tri(substitutions)],
                         substitutions[upper.tri(substitutions)]),
                 pred = c(sub_pred[lower.tri(sub_pred)],
                          sub_pred[upper.tri(sub_pred)]))

lims <- range(c(sub_df$obs, sub_df$pred))

sub_df %>%
    ggplot(aes(pred, obs)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "firebrick2") +
    scale_x_continuous("Predicted substitution count",
                       limits = lims) +
    scale_y_continuous("Observed substitution count",
                       limits = lims) +
    theme_classic()






# =================================================================*
# =================================================================*

# GAMMAS ----

# =================================================================*
# =================================================================*

ref <- with(pars, ref_genome$new(jackalope:::make_ref_genome(seqs[1:4])))
# Only substitutions plus gammas
mevo_ <- make_mevo(ref,
                   sub = list(model = "TN93", alpha_1 = pars$alpha_1,
                              alpha_2 = pars$alpha_2, beta = pars$beta,
                              pi_tcag = pars$pi_tcag),
                   site_var = pars$mats[1:ref$n_seqs()],
                   chunk_size = 100)
# Takes ~14 sec
var_set <- create_variants(ref, method = "phylo", method_info = tree,
                           mevo_obj = mevo_, n_cores = 4)

pos_df <- map_dfr(0:(ref$n_seqs()-1),
    function(s) {
        Z <- jackalope:::examine_mutations(var_set_ptr = var_set$genomes,
                                        var_ind = 0, seq_ind = s) %>%
            .[["pos"]]
        return(tibble(seq = as.integer(s), pos = as.integer(Z)))
    })

gamm_df <- map_dfr(0:(ref$n_seqs()-1),
    function(s) {
        gamma_mat <- mevo_$gamma_mats[[s+1]]
        Z <- pos_df[pos_df$seq == s, ][["pos"]]
        # Expected # mutations per gamma region:
        pred_mut = length(Z) / nrow(gamma_mat)
        df_ <- tibble(seq = s,
                      gamma = gamma_mat[,2],
                      pos = gamma_mat[,1],
                      obs = jackalope:::table_gammas(gamma_mat[,1], Z))
        df_$pred <- (df_$gamma / sum(df_$gamma)) * sum(df_$obs)
        return(df_)
    })


lims <- range(c(gamm_df$obs, gamm_df$pred))

gamm_df %>%
    ggplot(aes(pred, obs)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "firebrick2") +
    scale_x_continuous("Predicted mutation count",
                       limits = lims) +
    scale_y_continuous("Observed mutation count",
                       limits = lims) +
    theme_classic()






# =================================================================*
# =================================================================*

# PHYLOGENY ----

# =================================================================*
# =================================================================*

# Phylogenetic tree:
set.seed(89415648)
tree <- ape::rcoal(pars$n_vars)
tree$edge.length <- tree$edge.length * 0.1

ref <- with(pars, ref_genome$new(jackalope:::make_ref_genome(seqs[1:4])))
# Only substitutions again
mevo_ <- make_mevo(ref,
                   sub = list(model = "TN93", alpha_1 = pars$alpha_1,
                              alpha_2 = pars$alpha_2, beta = pars$beta,
                              pi_tcag = pars$pi_tcag))

var_set <- create_variants(ref, method = "phylo", method_info = tree,
                           mevo_obj = mevo_, n_cores = 4)


mutations <- lapply(0:(pars$n_vars-1),
                    function(v) {
                        Z <- jackalope:::view_mutations(var_set_ = var_set$genomes,
                                                     var_ind = v)
                        Zseq <- Z[,"seq"]
                        Z <- Z[, c("old_pos", "size_mod", "nucleos")]
                        ZZ <- apply(Z, 1, function(xx) gsub(" ", "",
                                                            paste(xx, collapse = "_")))
                        return(split(ZZ, Zseq))
                    })



shared_muts <- rep(list(matrix(0L, pars$n_vars, pars$n_vars)), ref$n_seqs())
for (s in 1:(ref$n_seqs())) {
    for (i in 1:(pars$n_vars-1)) {
        for (j in (i+1):(pars$n_vars)) {
            muts1 <- mutations[[i]][[s]]
            muts2 <- mutations[[j]][[s]]
            shared_muts[[s]][i,j] <- as.integer(sum(muts2 %in% muts1))
            shared_muts[[s]][j,i] <- as.integer(sum(muts1 %in% muts2))
        }
    }
    diag(shared_muts[[s]]) <- sapply(mutations, function(x) length(x[[s]]))
}; rm(s, i, j, muts1, muts2)


# "Mutational" distances between all pairs:
mut_dists <- rep(list(matrix(0L, pars$n_vars, pars$n_vars)), ref$n_seqs())
for (s in 1:(ref$n_seqs())) {
    for (i in 1:(pars$n_vars-1)) {
        for (j in (i+1):pars$n_vars) {
            # Total non-shared mutations for each:
            mut_dists[[s]][i,j] <- shared_muts[[s]][i,i] - shared_muts[[s]][i,j]
            mut_dists[[s]][j,i] <- shared_muts[[s]][j,j] - shared_muts[[s]][j,i]
        }
    }
}; rm(s, i, j)


# Phylogenetic distances between all pairs:
phylo_dists <- ape::cophenetic.phylo(tree)
rownames(phylo_dists) <- colnames(phylo_dists) <- NULL
# Dividing by two to get branch lengths between MRCA:
phylo_dists <- phylo_dists / 2


# Average mutation rate for ancestral sequences times total bp in each sequence:
mu_ <- mean(mevo_$q()) * mean(nchar(pars$seqs[1:4]))

# Expected mutational distances:
pred_mut_dists <- mu_ * phylo_dists


phy_df <- tibble(obs = lapply(mut_dists,
                              function(x) c(x[lower.tri(x)], x[upper.tri(x)])) %>%
                     do.call(what = c),
                 pred = rep(c(pred_mut_dists[lower.tri(pred_mut_dists)],
                              pred_mut_dists[upper.tri(pred_mut_dists)]),
                            ref$n_seqs()))

lims <- range(c(phy_df$obs, phy_df$pred))

phy_df %>%
    ggplot(aes(pred, obs)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "firebrick2") +
    scale_x_continuous("Predicted mutation count",
                       limits = lims) +
    scale_y_continuous("Observed mutation count",
                       limits = lims) +
    theme_classic()



# =================================================================*
# =================================================================*

# RATE METHOD ----

# =================================================================*
# =================================================================*

# Phylogenetic tree:
tree <- ape::rcoal(2)
tree$edge.length <- tree$edge.length * 0.1

ref <- with(pars, ref_genome$new(jackalope:::make_ref_genome(seqs[1:4])))

# Only substitutions
mevo_ <- make_mevo(ref,
                   sub = list(model = "TN93", alpha_1 = pars$alpha_1,
                              alpha_2 = pars$alpha_2, beta = pars$beta,
                              pi_tcag = pars$pi_tcag),
                   site_var = pars$mats[1:ref$n_seqs()])
var_set <- create_variants(ref, method = "phylo", method_info = tree,
                           mevo_obj = mevo_)

# Full sequences:
var_seqs <- lapply(0:(var_set$n_vars()-1),
                   function(v) jackalope:::view_var_genome(var_set_ = var_set$genomes,
                                                        var_ind = v))

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
    rate_cpp <- jackalope:::test_rate(start = start-1, end = end-1,
                                   var_ind = var_ind-1, seq_ind = seq_ind-1,
                                   var_set_ = var_set$genomes,
                                   sampler_ = jackalope:::mevo_obj_to_ptr(mevo_),
                                   gamma_mat_ = gamma_mat_new)
    rate_r <- get_seq_rate(var_seqs[[var_ind]][seq_ind], rates, start, end, gamma_mat_new)
    if (verbose & !isTRUE(all.equal(rate_cpp, rate_r))) {
        cat(sprintf("%i, %i, %i, %i\n", var_ind, seq_ind, start, end))
    }
    return(cbind(rate_cpp, rate_r))
}


mat <- matrix(0, var_set$n_vars() * ref$n_seqs(), 2)
i <- 1
for (var_ind in 1:var_set$n_vars()) {
    for (seq_ind in 1:ref$n_seqs()) {
        mat[i,] <- compare_rates(var_ind, seq_ind, rates = mevo_$q())
        i <- i + 1
    }
}


# Should be `TRUE`:
all.equal(mat[,1], mat[,2])

