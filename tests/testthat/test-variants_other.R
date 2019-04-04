

context("Testing basics of creating variants")

# library(jackalope)
# library(testthat)


reference <- create_genome(3, 100)
mevo_obj <- make_mevo(reference, list(model = "JC69", lambda = 0.1),
                      list(rate = 0.1, max_length = 10),
                      list(rate = 0.1, max_length = 10),
                      chunk_size = 0)



test_that("basics of variant creation work with theta method", {

    vars <- create_variants(reference, "theta", list(theta = 0.1, n_vars = 4), mevo_obj)

    expect_identical(vars$n_seqs(), reference$n_seqs())
    expect_identical(vars$n_vars(), 4L)

    vars2 <- create_variants(reference, "theta", list(theta = 4, n_vars = 4), mevo_obj)

    expect_gt(nrow(jackalope:::view_mutations(vars2$genomes, 0)),
              nrow(jackalope:::view_mutations(vars$genomes, 0)))

    expect_error(create_variants(reference, "theta", "theta", mevo_obj),
                 regexp = paste("argument `method_info` must be a named list or",
                                "numeric vector, with the names \"theta\" and",
                                "\"n_vars\". \"theta\" must be a single number > 0, and",
                                "\"n_vars\" must be a single whole number >= 1"))

})





test_that("basics of variant creation work with phylo method", {

    tr <- ape::rcoal(4)
    tr$edge.length <- tr$edge.length * 0.01

    vars <- create_variants(reference, "phylo", tr, mevo_obj)

    expect_identical(vars$n_seqs(), reference$n_seqs())
    expect_identical(vars$n_vars(), 4L)

    tr$edge.length <- tr$edge.length * 100

    vars2 <- create_variants(reference, "phylo", tr, mevo_obj)

    expect_gt(nrow(jackalope:::view_mutations(vars2$genomes, 0)),
              nrow(jackalope:::view_mutations(vars$genomes, 0)))

    expect_error(create_variants(reference, "phylo", "tr", mevo_obj),
                 regexp = paste("argument `method_info` must be of class \"phylo\",",
                                "\"multiPhylo\", or a list of \"phylo\" objects when",
                                "`method` = \"phylo\""))

})



test_that("basics of variant creation work with newick method", {

    tr <- ape::rcoal(4)

    tr_file <- paste0(tempdir(), "/test.tree")

    ape::write.tree(tr, tr_file)

    vars <- create_variants(reference, "newick", tr_file, mevo_obj)

    expect_identical(vars$n_seqs(), reference$n_seqs())
    expect_identical(vars$n_vars(), 4L)

    expect_error(create_variants(reference, "newick", tr, mevo_obj),
                 regexp = paste("argument `method_info` must be a single string",
                                "or a vector of strings of the same length as the",
                                "number of sequences, when `method` = \"newick\""))

})




test_that("errors occur when nonsense is input to create_variants", {

    expect_error(create_variants(reference, method = "coal_trees", method_info = c(),
                                 mevo_obj),
                 regexp = paste("argument `method_info` must be a single string or",
                                "a list when `method` = \"coal_trees\""))

    expect_error(create_variants(reference, method = "coal_sites", method_info = c(),
                                 mevo_obj),
                 regexp = paste("argument `method_info` must be a single string or",
                                "a list \\(for method = \"coal_sites\"\\)"))

})


