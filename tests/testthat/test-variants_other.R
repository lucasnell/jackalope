

# library(jackalope)
# library(testthat)

context("Testing basics of creating variants")


arg_list <- list(reference = create_genome(3, 100),
                 sub = sub_JC69(0.1, gamma_shape = 1))
arg_list$ins <- indels(rate = 0.1, max_length = 10)
arg_list$del <- indels(rate = 0.1, max_length = 10)


cv <- function(vars_info, al = arg_list) {
    arg_list_ <- c(list(vars_info = vars_info), al)
    vars <- do.call(create_variants, arg_list_)
    return(vars)
}

test_that("nonsense `sub` arg throws error", {
    al2 <- arg_list
    al2$sub <- "sub"
    expect_error(cv(vars_theta(0.1, n_vars = 4), al2),
                 regexp = "argument `sub` must be NULL or a \"sub_info\" object")
})



# vars_theta -----
test_that("basics of vars_theta work", {

    vi <- vars_theta(0.1, n_vars = 4)
    vars <- cv(vi)

    expect_identical(vars$n_chroms(), arg_list$reference$n_chroms())
    expect_identical(vars$n_vars(), 4L)

    vars2 <- cv(vars_theta(4, n_vars = 4))

    muts <- jackalope:::view_mutations(vars$ptr(), 0)
    muts2 <- jackalope:::view_mutations(vars2$ptr(), 0)

    expect_gt(sum(abs(muts2$size_mod)) + sum(muts2$size_mod == 0),
              sum(abs(muts$size_mod)) + sum(muts$size_mod == 0))

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

    expect_identical(vars$n_chroms(), arg_list$reference$n_chroms())
    expect_identical(vars$n_vars(), 4L)

    tr$edge.length <- tr$edge.length * 100

    vars2 <- cv(vars_phylo(tr))

    expect_gt(nrow(jackalope:::view_mutations(vars2$ptr(), 0)),
              nrow(jackalope:::view_mutations(vars$ptr(), 0)))

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

    expect_identical(vars$n_chroms(), arg_list$reference$n_chroms())
    expect_identical(vars$n_vars(), 4L)

    expect_error(vars_phylo(fn = tr),
                 regexp = "argument `fn` must be NULL or a character vector")

})




# basic output -----
test_that("basic diagnostic functions work for variants", {

    vars <- cv(vars_theta(theta = 0.1, n_vars = 4))

    Z <- jackalope:::examine_mutations(var_set_ptr = vars$ptr(),
                                       var_ind = 0, chrom_ind = 0)

    expect_identical(length(Z$pos), as.integer(sum(sapply(c("sub", "ins", "del"),
                                                          function(x) sum(Z[[x]])))))

})

