
# library(jackalope)
# library(testthat)

context("Testing R class methods and info")


n_chroms <- 10L
n_vars <- 5L
n_muts <- 100L
len <- 100L
len_sd <- 10.0

# Extract vector of chromosome strings:
chroms <- jackalope:::rando_chroms(n_chroms, len, len_sd, pi_tcag = c(8, 4, 2, 1))



# reference basics -----

test_that("Chromosomes from `create_genome` aren't very different from expectation.", {
    all_chroms <- paste(chroms, collapse = "")
    freq_obs <- sapply(c("T", "C", "A", "G"),
                       function(char) {
                           s2 <- gsub(char,"",all_chroms)
                           return((nchar(all_chroms) - nchar(s2)) / nchar(all_chroms))
                       })
    freq_obs <- as.numeric(freq_obs)
    freq_exp <- c(8, 4, 2, 1) / sum(c(8, 4, 2, 1))
    expect_identical(rank(freq_obs), rank(freq_exp))
})



test_that("initialization of ref_genome class with nonsense produces error", {
    expect_error(ref_genome$new("nonsense"),
                 paste("When initializing a ref_genome object, you need to use",
                       "an externalptr object"))
})

# Reference genome
ref <- ref_genome$new(jackalope:::make_ref_genome(chroms))
test_that("ref_genome class starts with the correct fields", {
    expect_is(ref$ptr(), "externalptr")
})


test_that("ref_genome print appears about right", {
    expect_output(print(ref),
                  paste0("< Set of ", n_chroms, " chromosomes >\n",
                         "\\# Total size: ", format(sum(nchar(chroms)),
                                                    big.mark = ","), " bp"))
})


test_that("ref_genome methods `set_names` and `clean_names` work properly", {

    og_names <- ref$chrom_names()
    unclean_names <- sprintf("nonsense names %i ';,\"", 1:ref$n_chroms())
    clean_names <- sprintf("nonsense_names_%i_____", 1:ref$n_chroms())

    ref$set_names(unclean_names)
    expect_identical(ref$chrom_names(), unclean_names)

    ref$clean_names()
    expect_identical(ref$chrom_names(), clean_names)

    ref$set_names(og_names)
    expect_identical(ref$chrom_names(), og_names)

})



test_that("ref_genome class methods produce correct output", {

    expect_identical(ref$n_chroms(), n_chroms)

    expect_identical(ref$sizes(), nchar(chroms))

    expect_identical(ref$chrom_names(), paste0("chrom", 1:length(chroms) - 1))

    expect_identical(jackalope:::view_ref_genome(ref$ptr()), chroms)

    nn <- paste0("__CHROM_",1:length(chroms))
    ref$set_names(nn)
    expect_identical(ref$chrom_names(), nn)

    ref$rm_chroms(nn[3])
    expect_identical(ref$chrom_names(), nn[-3])

    ref$add_chroms(chroms[3], nn[3])
    expect_identical(ref$chrom_names(), c(nn[-3], nn[3]))

    ref$merge_chroms()
    expect_identical(nchar(ref$chrom(1)), nchar(paste(c(chroms[-3], chroms[3]),
                                                         collapse = "")))

    nchars <- nchar(chroms)
    # Making sure of no removal when it shouldn't:
    ref <<- ref_genome$new(jackalope:::make_ref_genome(chroms))
    ref$filter_chroms(min(nchar(chroms)) - 1, "size")
    expect_identical(ref$n_chroms(), n_chroms, label = "after lack of size filtering")
    ref$filter_chroms(1 - 0.99 * min(nchars) / sum(nchars), "prop")
    expect_identical(ref$n_chroms(), n_chroms, label = "after lack of filtering")
    # Now seeing if it `filter_chroms` removes when it should
    ref$filter_chroms(min(nchars) + 1, "size")
    expect_lt(ref$n_chroms(), length(chroms))

    ref <- ref_genome$new(jackalope:::make_ref_genome(chroms))
    ref$filter_chroms(0.99 * sum(nchars[nchars > min(nchars)]) / sum(nchars), "prop")
    expect_lt(ref$n_chroms(), length(chroms))
})
# Restart object:
ref <- ref_genome$new(jackalope:::make_ref_genome(chroms))







# ============================================================`
# ============================================================`

# variants basics -----

# ============================================================`
# ============================================================`



# Create variants:
phy <- ape::rcoal(n_vars)
vars <- create_variants(ref, vars_phylo(phy), sub = sub_JC69(0.01))

test_that("variants class starts with the correct fields", {
    expect_is(vars$ptr(), "externalptr")
})


test_that("variants print appears about right", {
    expect_output(print(vars),
                  paste0("<< Variants object >>\n",
                         "\\# Variants: ", n_vars, "\n",
                         "\\# Mutations: "))
})



test_that("variants class methods", {

    expect_equal(vars$n_chroms(), length(chroms))

    expect_equal(vars$n_vars(), n_vars)

    # No indels, so these should be true:
    for (i in 1:n_vars) expect_equal(vars$sizes(i), nchar(chroms))

    expect_identical(vars$chrom_names(), paste0("chrom", 1:length(chroms)-1))

    expect_identical(vars$var_names(), phy$tip.label)

    # Variant chromosomes:
    var_chroms <- lapply(1:n_vars,
                       function(v) {
                           jackalope:::view_var_genome(vars$ptr(), v-1)
                       })
    vs0 <- lapply(1:n_vars, function(v) sapply(1:length(chroms),
                                               function(s) vars$chrom(v, s)))
    expect_identical(var_chroms, vs0)

    nv <- paste0("__VARS_", 1:n_vars)
    vars$set_names(nv)
    expect_identical(vars$var_names(), nv)

    vars$rm_vars(nv[1:2])
    expect_identical(vars$var_names(), nv[-1:-2])
})





# ============================================================`
# ============================================================`

# manual mutations -----

# ============================================================`
# ============================================================`


# Make empty var_set to compare mutations
vars <- variants$new(jackalope:::make_var_set(ref$ptr(), n_vars), ref$ptr())
vars_R <- replicate(n_vars, chroms, simplify = FALSE)


for (v in 1:n_vars) {
    ts <- vars_R[[v]]
    for (s in 1:length(chroms)) {
        m = 0;
        max_size = vars$sizes(v)[s]
        while (m < n_muts && max_size > 0) {
            pos = as.integer(runif(1) * max_size) + 1
            rnd = runif(1);
            if (rnd < 0.5) {
                str = jackalope:::rando_chroms(1, 1)
                if (nchar(str) != 1) stop("Improper size in sub")
                vars$add_sub(v, s, pos, str)
                substr(ts[s], pos, pos) <- str
            } else if (rnd < 0.75) {
                size = as.integer(rexp(1, 2.0) + 1.0)
                if (size > 10) size = 10
                str = jackalope:::rando_chroms(1, size)
                if (nchar(str) != size) stop("Improper size in insertion")
                vars$add_ins(v, s, pos, str)
                ts[s] <- paste0(substr(ts[s], 1, pos), str,
                                substr(ts[s], pos + 1, nchar(ts[s])))
                max_size <- max_size + size
            } else {
                size = as.integer(rexp(1, 2.0) + 1.0)
                if (size > 10) size = 10
                if ((pos + size - 1) > max_size) size = max_size - pos + 1;
                vars$add_del(v, s, pos, size)
                ts[s] <- paste0(substr(ts[s], 1, pos - 1),
                                substr(ts[s], pos + size, nchar(ts[s])))
                max_size <- max_size - size
            }
            m <- m + 1
            prev_rnd <- rnd
        }
    }
    vars_R[[v]] <- ts
}

# Converting to list like vars_R:
vars_cpp <- lapply(1:n_vars, function(v) sapply(1:n_chroms, function(s) vars$chrom(v, s)))

test_that("Mutations produced are accurate", {
    expect_identical(vars_R, vars_cpp)
})




# ============================================================`
# ============================================================`

# replace_Ns, gc/nt_prop -----

# ============================================================`
# ============================================================`


# Testing that replace_Ns works
chroms <- c("CCAANNNGG", "NNTTCCAAGG", "AACCTTGGGGGNNNNNN")
ref <- ref_genome$new(jackalope:::make_ref_genome(chroms))
ref$replace_Ns(c(1,0,0,0))

test_that("Replacing Ns works as predicted", {
    expect_identical(ref$chrom(1), "CCAATTTGG")
    expect_identical(ref$chrom(2), "TTTTCCAAGG")
    expect_identical(ref$chrom(3), "AACCTTGGGGGTTTTTT")
})


# Testing that gc_prop and nt_prob work
ref <- ref_genome$new(jackalope:::make_ref_genome(
    c(paste(c(rep("T", 50), rep("C", 50), rep("A", 50), rep("G", 50)), collapse = ""),
      paste(c(rep("T", 100), rep("C", 50), rep("A", 25), rep("G", 25)), collapse = ""),
      paste(c(rep("T", 50), rep("C", 25), rep("A", 25), rep("G", 100)), collapse = ""),
      paste(c(rep("T", 25), rep("C", 25), rep("A", 100), rep("G", 50)), collapse = ""))
    ))

test_that("gc_prop and nt_prob work for ref_genome as predicted", {

    expect_equal(ref$gc_prop(1, 1, 200), 100 / 200)
    expect_equal(ref$gc_prop(2, 1, 200), 75 / 200)
    expect_equal(ref$gc_prop(3, 1, 200), 125 / 200)
    expect_equal(ref$gc_prop(4, 1, 200), 75 / 200)

    expect_equal(ref$gc_prop(1, 1, 100), 50 / 100)
    expect_equal(ref$gc_prop(2, 1, 100), 0 / 100)
    expect_equal(ref$gc_prop(3, 1, 100), 25 / 100)
    expect_equal(ref$gc_prop(4, 1, 100), 25 / 100)


    expect_equal(ref$nt_prop('T', 1, 1, 200), 50 / 200)
    expect_equal(ref$nt_prop('T', 2, 1, 200), 100 / 200)
    expect_equal(ref$nt_prop('T', 3, 1, 200), 50 / 200)
    expect_equal(ref$nt_prop('T', 4, 1, 200), 25 / 200)

    expect_equal(ref$nt_prop('T', 1, 1, 100), 50 / 100)
    expect_equal(ref$nt_prop('T', 2, 1, 100), 100 / 100)
    expect_equal(ref$nt_prop('T', 3, 1, 100), 50 / 100)
    expect_equal(ref$nt_prop('T', 4, 1, 100), 25 / 100)

})




vars <- variants$new(jackalope:::make_var_set(ref$ptr(), 2), ref$ptr())

test_that("gc_prop and nt_prob work for variants as predicted", {

    expect_equal(vars$gc_prop(1, 1, 1, 200), 100 / 200)
    expect_equal(vars$gc_prop(1, 2, 1, 200), 75 / 200)
    expect_equal(vars$gc_prop(1, 3, 1, 200), 125 / 200)
    expect_equal(vars$gc_prop(1, 4, 1, 200), 75 / 200)

    expect_equal(vars$gc_prop(1, 1, 1, 100), 50 / 100)
    expect_equal(vars$gc_prop(1, 2, 1, 100), 0 / 100)
    expect_equal(vars$gc_prop(1, 3, 1, 100), 25 / 100)
    expect_equal(vars$gc_prop(1, 4, 1, 100), 25 / 100)


    expect_equal(vars$nt_prop('T', 1, 1, 1, 200), 50 / 200)
    expect_equal(vars$nt_prop('T', 1, 2, 1, 200), 100 / 200)
    expect_equal(vars$nt_prop('T', 1, 3, 1, 200), 50 / 200)
    expect_equal(vars$nt_prop('T', 1, 4, 1, 200), 25 / 200)

    expect_equal(vars$nt_prop('T', 1, 1, 1, 100), 50 / 100)
    expect_equal(vars$nt_prop('T', 1, 2, 1, 100), 100 / 100)
    expect_equal(vars$nt_prop('T', 1, 3, 1, 100), 50 / 100)
    expect_equal(vars$nt_prop('T', 1, 4, 1, 100), 25 / 100)

    # Adding some substitutions
    for (i in 1:100) vars$add_sub(2, 1, i, 'C')
    for (i in 101:200) vars$add_sub(2, 1, i, 'T')

    expect_equal(vars$gc_prop(2, 1, 1, 200), 100 / 200)
    expect_equal(vars$gc_prop(2, 1, 1, 100), 100 / 100)
    expect_equal(vars$gc_prop(2, 1, 101, 200), 0 / 100)

    expect_equal(vars$nt_prop('C', 2, 1, 1, 100), 100 / 100)
    expect_equal(vars$nt_prop('T', 2, 1, 101, 200), 100 / 100)

})







# ============================================================`
# ============================================================`

# +/- variants -----

# ============================================================`
# ============================================================`

ref <- ref_genome$new(jackalope:::make_ref_genome(chroms))

test_that("adding/removing/duplicating variants works as predicted", {

    vars0 <- create_variants(ref, NULL)
    expect_identical(vars0$n_vars(), 0L)

    vars0$add_vars("var0")
    expect_identical(vars0$n_vars(), 1L)
    expect_identical(sapply(1:ref$n_chroms(), function(i) vars0$chrom(1, i)),
                     sapply(1:ref$n_chroms(), function(i) ref$chrom(i)))

    vars1 <- create_variants(ref, vars_theta(0.1, 4), sub_JC69(0.001))
    vars1$dup_vars(vars1$var_names()[1:2])
    expect_identical(vars1$n_vars(), 6L)
    expect_identical(sapply(1:ref$n_chroms(), function(i) vars1$chrom(1, i)),
                     sapply(1:ref$n_chroms(), function(i) vars1$chrom(5, i)))
    expect_identical(sapply(1:ref$n_chroms(), function(i) vars1$chrom(2, i)),
                     sapply(1:ref$n_chroms(), function(i) vars1$chrom(6, i)))

})


