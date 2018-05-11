context("Testing molecular evolution accuracy")

library(gemino)

# Construct sequence with known number of each nucleotide
seq <- c(rep("T", 0.25e6), rep("C", 0.25e6), rep("A", 0.25e6), rep("G", 0.25e6))
seq <- sample(seq)
seq <- gemino:::cpp_merge_str(seq)

# Set molecular evolution parameters
N_ <- 1e3
pi_tcag_ = 0.1*1:4
alpha_1_ = 1
alpha_2_ = 2
beta_ = 0.5
xi_ = 0.25
psi_ = 1

Q <- matrix(beta_, 4, 4)
Q[2,1] <- Q[1,2] <- alpha_1_
Q[4,3] <- Q[3,4] <- alpha_2_
for (i in 1:4) Q[,i] <- Q[,i] * pi_tcag_[i]
diag(Q) <- 0

q <- rowSums(Q) + xi_

gm <- matrix(1L, 1e3, 2)
gm[,1] <- seq(1e3, 1e6, 1e3)
set.seed(1087437801)
gm[,2] <- rgamma(1e3, shape = 1, rate = 1)


# Make a matrix have ten columns
make_ten <- function(x) cbind(x, matrix(0, nrow(x), 10 - ncol(x)))


test_samp <- function() {
    vars <- gemino:::make_vars(seq, 1)
    gemino:::test_sampling(vs_sexp = vars,
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

muts <- lapply(1:nreps, function(i) test_samp())



# Checking the indel accuracy:

all_ins <- do.call(rbind, lapply(1:nreps, function(i) setNames(colSums(muts[[i]]$ins),
                                                               1:10)))
all_del <- do.call(rbind, lapply(1:nreps, function(i) setNames(colSums(muts[[i]]$del),
                                                               1:10)))
# Function to calculate the expected number of insertions or deletions of a given size:
indel_p <- function(ss) (N_ * 4 * xi_ * 0.5 / sum(q)) * (exp(-ss) / sum(exp(-1:-10)))


test_that("molecular evolution produces correct proportion of indel sizes", {
    expect_true(t.test(rowSums(all_ins) / rowSums(all_del), mu = 1)$p.value > 0.05)
    for (i in 1:10) {
        if (sum(all_ins[,i]) > 0) {
            p <- t.test(all_ins[,i], mu = indel_p(i))$p.value
        } else {
            p <- binom.test(0, N_ * nreps, p = indel_p(i) / N_)$p.value
        }
        expect_true(p > 0.05)
    }
    for (i in 1:10) {
        if (sum(all_del[,i]) > 0) {
            p <- t.test(all_del[,i], mu = indel_p(i))$p.value
        } else {
            p <- binom.test(0, N_ * nreps, p = indel_p(i) / N_)$p.value
        }
        expect_true(p > 0.05)
    }
})



# Checking the substition accuracy:


sub_df <- map_dfr(1:nreps, ~ data_frame(r = .x,
                                        t = muts[[.x]]$sub[1,],
                                        c = muts[[.x]]$sub[2,],
                                        a = muts[[.x]]$sub[3,],
                                        g = muts[[.x]]$sub[4,],
                                        to = c('t', 'c', 'a', 'g'))) %>%
    gather('from', 'count', t:g) %>%
    mutate(to = factor(to, levels = c('t', 'c', 'a', 'g')),
           from = factor(from, levels = c('t', 'c', 'a', 'g'),
                         labels = paste("from:", c('t', 'c', 'a', 'g'))))

sub_df <- sub_df %>%
    filter(!((from == "from: t" & to == "t") |
                 (from == "from: c" & to == "c") |
                 (from == "from: a" & to == "a") |
                 (from == "from: g" & to == "g")))

sub_df %>%
    ggplot(aes(to, count, color = to)) +
    geom_point(position = position_jitter(width = 0.25, height = 0), shape = 1) +
    theme_classic() +
    facet_wrap(~ from, nrow = 2, scales = "free") +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    geom_point(data = sub_df %>%
                   distinct(to, from) %>%
                   arrange(to, from) %>%
                   mutate(count = map2_dbl(from, to,
                                           function(ii,jj) {
                                               i <- as.integer(ii)
                                               j <- as.integer(jj)
                                               e <- N_ * ((q[i] - xi_) / sum(q)) *
                                                   (Q[i, j] / sum(Q[i,]))
                                               return(e)
                                           })),
               color = 'black', shape = 19, size = 3)



# Now looking at positions

pos_df <- map_dfr(1:nreps, function(i) {
    data_frame(rep = i,
               gamma = gm[,2],
               pos = gm[,1],
               count = gemino:::table_gammas(gm[,1], muts[[i]]$pos))
})


