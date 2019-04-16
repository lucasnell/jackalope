
#'
#' The files `files/ms_out*.txt` are used for these tests, and all but
#' `ms_out.txt` have errors intentionally added:
#'
#' - `ms_out_err1.txt`: no sites and no trees
#' - `ms_out_err2.txt`: a position is removed in the seg. sites
#'   (`"error bad site pos"`),
#'   a variant is removed in a tree (`"error missing tree tip"`).
#' - `ms_out_err3.txt`: a position is removed in a tree (`"error missing tree pos"`),
#'   a site is removed in seg. sites
#'
#'
#'



context("Testing creating variants from coalescent simulations")

# library(jackalope)
# library(testthat)

# necessary objects -----

arg_list <- list(reference = create_genome(3, 100), sub = sub_JC69(0.1),
                 ins = indels(rate = 0.1, max_length = 10),
                 del = indels(rate = 0.1, max_length = 10))






# ==============================================================================`
# ==============================================================================`

# MS OBJECTS -----

# ==============================================================================`
# ==============================================================================`



#' Check for proper output and errors when creating variants from coalescent objects
#'
#' coal_obj must be a coalescent object from scrm or coala.


coal_obj_run_trees <- function(coal_obj, pkg) {

    arg_list_ <- c(list(vars_info = vars_gtrees(obj = coal_obj)), arg_list)
    vars <- do.call(create_variants, arg_list_)
    expect_equal(vars$n_vars(), 5, info = paste("works with gene trees -", pkg))

    invisible(NULL)
}

coal_obj_run_sites <- function(coal_obj, pkg) {
    arg_list_ <- c(list(vars_info = vars_ssites(obj = coal_obj)), arg_list)
    vars <- do.call(create_variants, arg_list_)
    expect_equal(vars$n_vars(), 5, info = paste("works with seg. sites -", pkg))
    invisible(NULL)
}
coal_obj_err_no_trees <- function(coal_obj, pkg) {
    coal_obj$trees <- NULL
    expect_error(
        vars_gtrees(obj = coal_obj),
        regexp = paste("`obj` must be \\(1\\) a list with a",
                       "`trees` field present or \\(2\\) a list of",
                       "lists, each sub-list containing a `trees`",
                       "field of length 1."),
        info = paste("returns error when no gene trees provided -", pkg))
    invisible(NULL)
}
coal_obj_err_no_sites <- function(coal_obj, pkg) {
    coal_obj$seg_sites <- NULL
    expect_error(
        vars_ssites(obj = coal_obj),
        regexp = paste("argument `obj` must be NULL or a list with a `seg_sites`",
                       "field present."),
        info = paste("returns error when no seg. sites provided -", pkg))
    invisible(NULL)
}
coal_obj_err_missing_tree_pos <- function(coal_obj, pkg) {
    i <- which(sapply(coal_obj$trees, length) > 1)[1]
    coal_obj$trees[[i]][1] <- substr(coal_obj$trees[[i]][1],
                                     regexpr("\\]", coal_obj$trees[[i]][1])[[1]]+1,
                                     9999)

    arg_list_ <- c(list(vars_info = vars_gtrees(obj = coal_obj)), arg_list)

    expect_error(
        do.call(create_variants, arg_list_),
        regexp = paste("A coalescent string appears to include",
                       "recombination but does not include sizes for each region."),
        info = paste("returns error when improper gene trees provided -", pkg))
    invisible(NULL)
}
coal_obj_err_miss_tree_tip <- function(coal_obj, pkg) {
    i <- which(sapply(coal_obj$trees, length) > 1)[1]
    strs <- coal_obj$trees[[i]]
    pos <- as.integer(sapply(strs, function(x) substr(x, 2, regexpr("\\]", x)[[1]]-1)))
    tr <- ape::read.tree(text = strs)
    tr <- lapply(tr, ape::drop.tip, tip = "1")
    class(tr) <- "multiPhylo"
    coal_obj$trees[[i]] <- sprintf("[%i]%s", pos, ape::write.tree(tr))

    arg_list_ <- c(list(vars_info = vars_gtrees(obj = coal_obj)), arg_list)

    expect_error(
        do.call(create_variants, arg_list_),
        regexp = "all gene trees must have the same number of tips.",
        info = paste("returns error when improper gene trees provided -", pkg))
    invisible(NULL)
}
coal_obj_err_bad_site_pos <- function(coal_obj, pkg) {
    i <- which(sapply(coal_obj$seg_sites, ncol) > 1)[1]
    mat <- coal_obj$seg_sites[[i]]
    mat <- as.matrix(mat)
    colnames(mat) <- c("XX", colnames(mat)[-1])
    coal_obj$seg_sites[[i]] <- mat
    expect_error(
        suppressWarnings({
            vars_ssites(obj = coal_obj)
        }),
        regexp = paste("Positions in one or more segregating-sites matrices",
                       "are producing NAs"),
        info = paste("returns error with improper seg. sites col names -", pkg))
    invisible(NULL)
}






# scrm checks -----

test_that("variant creation works with scrm coalescent object", {

    skip_if_not_installed("scrm")
    library(scrm)
    coal_obj <- scrm("5 3 -r 3.1 100 -t 10 -T -L")

    coal_obj_run_trees(coal_obj, "scrm")
    coal_obj_run_sites(coal_obj, "scrm")
    coal_obj_err_no_trees(coal_obj, "scrm")
    coal_obj_err_no_sites(coal_obj, "scrm")
    coal_obj_err_missing_tree_pos(coal_obj, "scrm")
    coal_obj_err_miss_tree_tip(coal_obj, "scrm")
    coal_obj_err_bad_site_pos(coal_obj, "scrm")


})


# coala checks -----

test_that("variant creation works with coala coalescent object", {

    skip_if_not_installed("coala")
    library(coala)
    model <- coal_model(sample_size = 5, loci_number = 3, loci_length = 100) +
        feat_recombination(2) +
        feat_mutation(5) +
        sumstat_trees() +
        sumstat_seg_sites()
    coal_obj <- simulate(model)

    coal_obj_run_trees(coal_obj, "coala")
    coal_obj_run_sites(coal_obj, "coala")
    coal_obj_err_no_trees(coal_obj, "coala")
    coal_obj_err_no_sites(coal_obj, "coala")
    coal_obj_err_missing_tree_pos(coal_obj, "coala")
    coal_obj_err_miss_tree_tip(coal_obj, "coala")
    coal_obj_err_bad_site_pos(coal_obj, "coala")

})








# ==============================================================================`
# ==============================================================================`

# MS FILES -----

# ==============================================================================`
# ==============================================================================`



test_that("variant creation works with ms-style file output", {


    .p <- function(x) test_path(sprintf("files/%s.txt", x))

    expect_equal({
        arg_list_ <- c(list(vars_info = vars_gtrees(fn = .p("ms_out"))), arg_list)
        vars <- do.call(create_variants, arg_list_)
        vars$n_vars()
    }, 5, info = paste("works with gene trees - ms-file"))

    expect_equal({
        arg_list_ <- c(list(vars_info = vars_gtrees(fn = .p("ms_out"))), arg_list)
        vars <- do.call(create_variants, arg_list_)
        vars$n_vars()
    }, 5, info = "works with seg. sites - ms-file")

    expect_error(
        vars_gtrees(fn = .p("ms_out_err1")),
        regexp = "one or more sequences have no trees",
        info = "returns error when no gene trees provided - ms-file")

    expect_error(
        vars_ssites(fn = .p("ms_out_err1")),
        regexp = paste("One or more seg. sites matrices from a ms-style",
                       "file output have no variant information specified."),
        info = "returns error when no seg. sites provided - ms-file")


    expect_error({
        arg_list_ <- c(list(vars_info = vars_gtrees(fn = .p("ms_out_err3"))), arg_list)
        do.call(create_variants, arg_list_)
    },
    regexp = paste("A coalescent string appears to include",
                   "recombination but does not include sizes for each region."),
    info = "returns error when improper gene trees provided - ms-file")


    expect_error({
        arg_list_ <- c(list(vars_info = vars_gtrees(fn = .p("ms_out_err3"))), arg_list)
        do.call(create_variants, arg_list_)
    },
    regexp = paste("A coalescent string appears to include",
                   "recombination but does not include sizes for each region."),
    info = "returns error when improper gene trees provided - ms-file")


    expect_error({
        arg_list_ <- c(list(vars_info = vars_gtrees(fn = .p("ms_out_err2"))), arg_list)
        do.call(create_variants, arg_list_)
    },
    regexp = "all gene trees must have the same number of tips.",
    info = "returns error when improper gene trees provided - ms-file")


    expect_error(
        vars_ssites(fn = test_path(paste0("files/", "ms_out_err2.txt"))),
        regexp = paste("the listed positions for each site \\(line starting",
                       "with 'positions:'\\) does not have a length that's",
                       "the same as the \\# sites as given by the line",
                       "starting with 'segsites:'"),
        info = paste("returns error with improper seg. sites positions - ms-file"))


    expect_error(
        vars_ssites(fn = test_path("files/ms_out_err3.txt")),
        regexp = paste("the listed number of sites \\(line starting with",
                       "'segsites:'\\) does not agree with the number of",
                       "items in the 3th line of segregating sites",
                       "info \\(ones filled with 0s and 1s\\)"),
        info = paste("returns error with improper seg. sites matrices - ms file"))


    reference2 <- create_genome(3, 1000)

    expect_error({
        arg_list_ <- c(list(vars_info = vars_gtrees(fn = .p("ms_out"))), arg_list)
        arg_list_$reference <- reference2
        do.call(create_variants, arg_list_)
    },
    regexp = paste("A coalescent string appears to include recombination",
                   "but the combined sizes of all regions don't match the",
                   "size of the sequence"))

})



test_that("variant creation returns error with improper ref_genome input", {
    .p <- function(x) test_path(sprintf("files/%s.txt", x))
    expect_error({
        arg_list_ <- c(list(vars_info = vars_gtrees(fn = .p("ms_out"))), arg_list)
        arg_list_$reference <- list(1, 2)
        do.call(create_variants, arg_list_)
    },
    regexp = paste("For the `create_variants` function in jackalope,",
                   "argument `reference` must be a \"ref_genome\" object."))
})

test_that("seg. sites produces the correct number of mutations", {
    ms_file <- test_path("files/ms_out.txt")
    # Have to make bigger ref. genome to make sure site positions don't overlap.
    reference2 <- create_genome(3, 100e3)
    arg_list_ <- c(list(vars_info = vars_ssites(fn = ms_file)), arg_list)
    arg_list_$reference <- reference2
    vars <- do.call(create_variants, arg_list_)

    msf <- readLines(ms_file)[-1:-2]
    msf <- msf[grepl("^0|^1", msf)]
    msf <- do.call(c, strsplit(msf, ""))
    n_muts_by_var <-
        sapply(0:4, function(i) nrow(jackalope:::view_mutations(vars$genomes, i)))
    expect_equal(sum(as.integer(msf)), sum(n_muts_by_var))
})





# ==============================================================================`
# ==============================================================================`

# GENERAL NONSENSE -----

# ==============================================================================`
# ==============================================================================`


test_that("errors occur when nonsense is input to vars_gtrees", {

    expect_error(vars_gtrees(NULL, NULL),
                 regexp = "either argument `obj` or `fn` must be provided")
    expect_error(vars_gtrees(1, 1),
                 regexp = "only one argument \\(`obj` or `fn`\\) should be provided")
    expect_error(vars_gtrees(1, NULL),
                 regexp = paste("argument `obj` must be \\(1\\) a list with a",
                                "`trees` field present or \\(2\\) a list of lists,",
                                "each sub-list containing a `trees` field of length 1"))
    expect_error(vars_gtrees(NULL, 1),
                 regexp = "argument `fn` must be NULL or a single string")

})
test_that("errors occur when nonsense is input to vars_ssites", {

    expect_error(vars_ssites(NULL, NULL),
                 regexp = "either argument `obj` or `fn` must be provided")
    expect_error(vars_ssites(1, 1),
                 regexp = "only one argument \\(`obj` or `fn`\\) should be provided")
    expect_error(vars_ssites(1, NULL),
                 regexp = paste("argument `obj` must be NULL or a list with a",
                                "`seg_sites` field present"))
    expect_error(vars_ssites(NULL, 1),
                 regexp = "argument `fn` must be NULL or a single string")

})
