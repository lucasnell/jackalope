

context("Testing basics of creating variants")

# library(jackalope)
# library(testthat)


arg_list <- list(reference = create_genome(3, 100),
                 sub = sub_JC69(0.1))
arg_list$ins <- indels(rate = 0.1, max_length = 10)
arg_list$del <- indels(rate = 0.1, max_length = 10)
arg_list$gamma_mats <- site_var(arg_list$reference, shape = 2, region_size = 10)


cv <- function(vars_info, al = arg_list) {
    arg_list_ <- c(list(vars_info = vars_info), al)
    vars <- do.call(create_variants, arg_list_)
    return(vars)
}

test_that("missing `sub` arg throws error", {
    al2 <- arg_list
    al2$sub <- NULL
    expect_error(cv(vars_theta(0.1, n_vars = 4), al2),
                 regexp = paste("argument `sub` must be provided if you want to create",
                                "variants using any method other than a VCF file"))
})


# vars_theta -----
test_that("basics of vars_theta work", {

    vars <- cv(vars_theta(0.1, n_vars = 4))

    expect_identical(vars$n_seqs(), arg_list$reference$n_seqs())
    expect_identical(vars$n_vars(), 4L)

    vars2 <- cv(vars_theta(4, n_vars = 4))

    expect_gt(nrow(jackalope:::view_mutations(vars2$genomes, 0)),
              nrow(jackalope:::view_mutations(vars$genomes, 0)))

    expect_error(vars_theta("0.1", n_vars = 4),
                 regexp = "argument `theta` must be a single number >= 0.")
    expect_error(vars_theta(0.1, n_vars = 1),
                 regexp = "argument `n_vars` must be a single integer >= 2.")
})



# vars_phylo w obj -----
test_that("basics of vars_phylo with object work", {

    tr <- ape::rcoal(4)
    tr$edge.length <- tr$edge.length * 0.01

    vars <- cv(vars_phylo(tr))

    expect_identical(vars$n_seqs(), arg_list$reference$n_seqs())
    expect_identical(vars$n_vars(), 4L)

    tr$edge.length <- tr$edge.length * 100

    vars2 <- cv(vars_phylo(tr))

    expect_gt(nrow(jackalope:::view_mutations(vars2$genomes, 0)),
              nrow(jackalope:::view_mutations(vars$genomes, 0)))

    expect_error(vars_phylo("tr"),
                 regexp = paste("argument `obj` must be NULL or of class \"phylo\",",
                       "\"multiPhylo\", or a list of \"phylo\" objects"))

})


# vars_phylo w file -----
test_that("basics of vars_phylo with file work", {

    tr <- ape::rcoal(4)

    tr_file <- paste0(tempdir(check = TRUE), "/test.tree")

    ape::write.tree(tr, tr_file)

    vars <- cv(vars_phylo(fn = tr_file))

    expect_identical(vars$n_seqs(), arg_list$reference$n_seqs())
    expect_identical(vars$n_vars(), 4L)

    expect_error(vars_phylo(fn = tr),
                 regexp = "argument `fn` must be NULL or a character vector")

})




# basic output -----
test_that("basic diagnostic functions work for variants", {

    vars <- cv(vars_theta(theta = 0.1, n_vars = 4))

    Z <- jackalope:::examine_mutations(var_set_ptr = vars$genomes,
                                       var_ind = 0, seq_ind = 0)

    expect_identical(length(Z$pos), as.integer(sum(sapply(c("sub", "ins", "del"),
                                                          function(x) sum(Z[[x]])))))
    expect_identical(jackalope:::table_gammas(seq(9, 99, 10), 0:99), rep(10, 10))

})

