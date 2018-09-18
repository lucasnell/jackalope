
context("Testing R class methods and info")

# library(gemino)
# library(testthat)


n_seqs <- 10
n_vars <- 5
n_muts <- 100
len <- 100
len_sd <- 10

# Extract vector of sequence strings:
seqs <- gemino:::rando_seqs(n_seqs, len, len_sd)


test_that(paste("Random sequences using `create_genome` aren't significantly different",
                "from expectation."),
          {
              # Frequency of `char` in string `s`
              char_frequencies <- function(char, s) {
                  s2 <- gsub(char,"",s)
                  return((nchar(s) - nchar(s2)) / nchar(s))
              }
              freq_obs <- t(sapply(1:n_seqs,
                                   function(i) sapply(c("T", "C", "A", "G"),
                                                      char_frequencies,
                                                      s = seqs[i])))
              freq_exp <- matrix(rep(0.25, n_seqs * 4), n_seqs)
              p <- chisq.test(freq_obs, freq_exp, simulate.p.value = TRUE)$p.value
              expect_gt(p, 0.05)
          }
)



# Reference genome
ref <- ref_genome$new(gemino:::make_ref_genome(seqs))
test_that("ref_genome class starts with the correct fields", {
    expect_is(ref$genome, "externalptr")
    expect_null(ref$digests)
})

test_that("ref_genome class methods", {

    expect_equal(ref$n_seqs(), n_seqs)

    expect_equal(ref$sizes(), nchar(seqs))

    expect_equal(ref$names(), paste0("seq", 1:length(seqs)-1))

    for (i in 1:n_seqs) expect_equal(ref$extract_seq(i), seqs[i])

    nn <- paste0("__SEQ_",1:length(seqs))
    ref$set_names(nn)
    expect_equal(ref$names(), nn)

    ref$rm_seqs(nn[3])
    expect_equal(ref$names(), nn[-3])

    ref$merge_seqs()
    expect_equal(nchar(ref$extract_seq(1)), nchar(paste(seqs[-3], collapse = "")))

    nchars <- nchar(seqs)
    # Making sure of no removal when it shouldn't:
    ref <<- ref_genome$new(gemino:::make_ref_genome(seqs))
    ref$filter_seqs(min(nchar(seqs)) - 1, "size")
    expect_equal(ref$n_seqs(), n_seqs, label = "after lack of size filtering")
    ref$filter_seqs(1 - 0.99 * min(nchars) / sum(nchars), "prop")
    expect_equal(ref$n_seqs(), n_seqs, label = "after lack of filtering")
    # Now seeing if it `filter_seqs` removes when it should
    ref$filter_seqs(min(nchars) + 1, "size")
    expect_lt(ref$n_seqs(), length(seqs))

    ref <<- ref_genome$new(gemino:::make_ref_genome(seqs))
    ref$filter_seqs(0.99 * sum(nchars[nchars > min(nchars)]) / sum(nchars), "prop")
    expect_lt(ref$n_seqs(), length(seqs))
})
# Restart object:
ref <- ref_genome$new(gemino:::make_ref_genome(seqs))

# Molecular evolution info:
mev <- make_mevo(ref, list(model = "JC69", lambda = 0.05))
# <test-make_mevo.R already tested this class>

# Create variants:
phy <- ape::rcoal(n_vars)
vars <- create_variants(ref, "phy", phy, mevo_obj = mev)
test_that("variants class starts with the correct fields", {
    expect_is(vars$genomes, "externalptr")
    expect_null(vars$digests)
})


test_that("variants class methods", {

    expect_equal(vars$n_seqs(), length(seqs))

    expect_equal(vars$n_vars(), n_vars)

    # No indels, so these should be true:
    for (i in 1:n_vars) expect_equal(vars$sizes(i), nchar(seqs))

    expect_equal(vars$seq_names(), paste0("seq", 1:length(seqs)-1))

    expect_equal(vars$var_names(), phy$tip.label)

    # Variant sequences:
    var_seqs <- lapply(1:n_vars,
                       function(v) {
                           gemino:::view_var_genome(vars$genomes, v-1)
                       })
    vs0 <- lapply(1:n_vars, function(v) sapply(1:length(seqs),
                                               function(s) vars$extract_seq(v, s)))
    expect_identical(var_seqs, vs0)

    nv <- paste0("__VARS_", 1:n_vars)
    vars$set_names(nv)
    expect_equal(vars$var_names(), nv)

    vars$rm_vars(nv[1:2])
    expect_equal(vars$var_names(), nv[-1:-2])
})



