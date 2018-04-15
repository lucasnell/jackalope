context("Testing mutation accuracy")

# Simulating mutations in R and gemino's C++ code and making sure they produce
# the same output

library(gemino)

n_seqs <- 10
n_vars <- 10
n_muts <- 500
len <- 100


seqs <- rando_seqs(n_seqs, len)

test_that("Random sequences using `rando_seqs` have the correct lengths.", {
    expect_equal(nchar(seqs), rep(len, n_seqs))
})

vars <- gemino:::make_vars(seqs, n_vars)
vars_R <- replicate(n_vars, seqs, simplify = FALSE)


for (v in 0:(n_vars-1)) {
    ts <- vars_R[[(v+1)]]
    for (s in 0:(length(seqs)-1)) {
        m = 0;
        max_size = gemino:::see_sizes(vars, v)[(s+1)]
        while (m < n_muts && max_size > 0) {
            pos = as.integer(runif(1) * max_size)
            rnd = runif(1);
            if (rnd < 0.5) {
                str = gemino:::rando_seqs(1, 1)
                if (nchar(str) != 1) stop("Improper size in sub")
                gemino:::add_substitution(vars, var_ind = v, seq_ind = s,
                                          nucleo = str, new_pos_ = pos)
                substr(ts[(s+1)], pos+1, pos+1) <- str
            } else if (rnd < 0.75) {
                size = as.integer(rexp(1, 2.0) + 1.0)
                if (size > 10) size = 10
                str = gemino:::rando_seqs(1, size)
                if (nchar(str) != size) stop("Improper size in insertion")
                gemino:::add_insertion(vars, var_ind = v, seq_ind = s,
                                       nucleos_ = str, new_pos_ = pos)
                ts[(s+1)] <- paste0(substr(ts[(s+1)], 1, pos + 1), str,
                                    substr(ts[(s+1)], pos + 2, nchar(ts[(s+1)])))
                max_size <- max_size + size
            } else {
                size = as.integer(rexp(1, 2.0) + 1.0)
                if (size > 10) size = 10
                if (size > (max_size - pos)) size = max_size - pos;
                gemino:::add_deletion(vars, var_ind = v, seq_ind = s,
                                      size_ = size, new_pos_ = pos)
                ts[(s+1)] <- paste0(substr(ts[(s+1)], 1, pos),
                                    substr(ts[(s+1)], pos + 1 + size,
                                           nchar(ts[(s+1)])))
                max_size <- max_size - size
            }
            m <- m + 1
            prev_rnd <- rnd
        }
    }
    vars_R[[(v+1)]] <- ts
}

# Converting to list like vars_R:
vars_cpp <- lapply(1:n_vars, function(i) gemino:::see_vg(vars, i-1))

# identical(vars_R, vars_cpp)

test_that("Mutations produced are accurate", {
    expect_identical(vars_R, vars_cpp)
})
