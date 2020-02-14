

# library(jackalope)
# library(testthat)

context("Testing basics of creating haplotypes")


arg_list <- list(reference = create_genome(3, 100),
                 sub = sub_JC69(0.1, gamma_shape = 1))
arg_list$ins <- indels(rate = 0.1, max_length = 10)
arg_list$del <- indels(rate = 0.1, max_length = 10)


cv <- function(haps_info, al = arg_list) {
    arg_list_ <- c(list(haps_info = haps_info), al)
    haps <- do.call(create_haplotypes, arg_list_)
    return(haps)
}

test_that("nonsense `sub` arg throws error", {
    al2 <- arg_list
    al2$sub <- "sub"
    expect_error(cv(haps_theta(0.1, n_haps = 4), al2),
                 regexp = "argument `sub` must be NULL or a \"sub_info\" object")
})



# haps_theta -----
test_that("basics of haps_theta work", {

    vi <- haps_theta(0.1, n_haps = 4)
    haps <- cv(vi)

    expect_identical(haps$n_chroms(), arg_list$reference$n_chroms())
    expect_identical(haps$n_haps(), 4L)

    haps2 <- cv(haps_theta(4, n_haps = 4))

    muts <- jackalope:::view_mutations(haps$ptr(), 0)
    muts2 <- jackalope:::view_mutations(haps2$ptr(), 0)

    expect_gt(sum(abs(muts2$size_mod)) + sum(muts2$size_mod == 0),
              sum(abs(muts$size_mod)) + sum(muts$size_mod == 0))

    expect_error(haps_theta("0.1", n_haps = 4),
                 regexp = "argument `theta` must be a single number >= 0.")
    expect_error(haps_theta(0.1, n_haps = 1),
                 regexp = "argument `n_haps` must be a single integer >= 2.")
})


# haps_theta -----
test_that("basics of haps_theta work - with exact indel simulation", {

    vi <- haps_theta(0.1, n_haps = 4)
    haps <- cv(vi, c(list(epsilon = 0), arg_list))

    expect_identical(haps$n_chroms(), arg_list$reference$n_chroms())
    expect_identical(haps$n_haps(), 4L)

    haps2 <- cv(haps_theta(4, n_haps = 4), c(list(epsilon = 0), arg_list))

    muts <- jackalope:::view_mutations(haps$ptr(), 0)
    muts2 <- jackalope:::view_mutations(haps2$ptr(), 0)

    expect_gt(sum(abs(muts2$size_mod)) + sum(muts2$size_mod == 0),
              sum(abs(muts$size_mod)) + sum(muts$size_mod == 0))

    expect_error(haps_theta("0.1", n_haps = 4),
                 regexp = "argument `theta` must be a single number >= 0.")
    expect_error(haps_theta(0.1, n_haps = 1),
                 regexp = "argument `n_haps` must be a single integer >= 2.")
})



# haps_phylo w obj -----
test_that("basics of haps_phylo with object work", {

    tr <- ape::rcoal(4)
    tr$edge.length <- tr$edge.length * 0.01

    haps <- cv(haps_phylo(tr))

    expect_identical(haps$n_chroms(), arg_list$reference$n_chroms())
    expect_identical(haps$n_haps(), 4L)

    tr$edge.length <- tr$edge.length * 100

    haps2 <- cv(haps_phylo(tr))

    expect_gt(nrow(jackalope:::view_mutations(haps2$ptr(), 0)),
              nrow(jackalope:::view_mutations(haps$ptr(), 0)))

    expect_error(haps_phylo("tr"),
                 regexp = paste("argument `obj` must be NULL or of class \"phylo\",",
                       "\"multiPhylo\", or a list of \"phylo\" objects"))

})


# haps_phylo w file -----
test_that("basics of haps_phylo with file work", {

    tr <- ape::rcoal(4)

    tr_file <- paste0(tempdir(check = TRUE), "/test.tree")

    ape::write.tree(tr, tr_file)

    haps <- cv(haps_phylo(fn = tr_file))

    expect_identical(haps$n_chroms(), arg_list$reference$n_chroms())
    expect_identical(haps$n_haps(), 4L)

    expect_error(haps_phylo(fn = tr),
                 regexp = "argument `fn` must be NULL or a character vector")

})




# basic output -----
test_that("basic diagnostic functions work for haplotypes", {

    haps <- cv(haps_theta(theta = 0.1, n_haps = 4))

    Z <- jackalope:::examine_mutations(hap_set_ptr = haps$ptr(),
                                       hap_ind = 0, chrom_ind = 0)

    expect_identical(length(Z$pos), as.integer(sum(sapply(c("sub", "ins", "del"),
                                                          function(x) sum(Z[[x]])))))

})

