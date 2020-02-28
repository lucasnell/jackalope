
#'
#' The files `files/ms_out*.txt` are used for these tests, and all but
#' `ms_out.txt` have errors intentionally added:
#'
#' - `ms_out_err1.txt`: no sites and no trees
#' - `ms_out_err2.txt`: a position is removed in the seg. sites
#'   (`"error bad site pos"`),
#'   a haplotype is removed in a tree (`"error missing tree tip"`).
#' - `ms_out_err3.txt`: a position is removed in a tree (`"error missing tree pos"`),
#'   a site is removed in seg. sites
#'
#'
#'




# library(jackalope)
# library(testthat)

context("Testing creating haplotypes from coalescent simulations")

# necessary objects -----

arg_list <- list(reference = create_genome(3, 100), sub = sub_JC69(0.1),
                 ins = indels(rate = 0.1, max_length = 10),
                 del = indels(rate = 0.1, max_length = 10))






# ==============================================================================`
# ==============================================================================`

# MS OBJECTS -----

# ==============================================================================`
# ==============================================================================`



#' Check for proper output and errors when creating haplotypes from coalescent objects
#'
#' coal_obj must be a coalescent object from scrm or coala.
#'

coal_obj_run_trees <- function(coal_obj, pkg) {

    arg_list_ <- c(list(haps_info = haps_gtrees(obj = coal_obj)), arg_list)
    haps <- do.call(create_haplotypes, arg_list_)
    expect_equal(haps$n_haps(), 5, info = paste("works with gene trees -", pkg))

    invisible(NULL)
}

coal_obj_run_sites <- function(coal_obj, pkg) {
    arg_list_ <- c(list(haps_info = haps_ssites(obj = coal_obj)), arg_list)
    haps <- do.call(create_haplotypes, arg_list_)
    expect_equal(haps$n_haps(), 5, info = paste("works with seg. sites -", pkg))
    invisible(NULL)
}
coal_obj_err_no_trees <- function(coal_obj, pkg) {
    coal_obj$trees <- NULL
    expect_error(
        haps_gtrees(obj = coal_obj),
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
        haps_ssites(obj = coal_obj),
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

    arg_list_ <- c(list(haps_info = haps_gtrees(obj = coal_obj)), arg_list)

    expect_error(
        do.call(create_haplotypes, arg_list_),
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

    arg_list_ <- c(list(haps_info = haps_gtrees(obj = coal_obj)), arg_list)

    expect_error(
        do.call(create_haplotypes, arg_list_),
        regexp = "One or more trees have differing tip labels.",
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
            haps_ssites(obj = coal_obj)
        }),
        regexp = paste("Positions in one or more segregating-sites matrices",
                       "are producing NAs"),
        info = paste("returns error with improper seg. sites col names -", pkg))
    invisible(NULL)
}






# scrm checks -----

test_that("haplotype creation works with scrm coalescent object", {

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


# scrm checks - exact -----

test_that("haplotype creation works with scrm coalescent object", {

    .arg_list <- arg_list

    arg_list <- c(list(epsilon = 0), arg_list)

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

    arg_list <- .arg_list

})


# coala checks -----

test_that("haplotype creation works with coala coalescent object", {

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



test_that("haplotype creation works with ms-style file output", {


    .p <- function(x) test_path(sprintf("files/%s.txt", x))

    expect_equal({
        arg_list_ <- c(list(haps_info = haps_gtrees(fn = .p("ms_out"))), arg_list)
        haps <- do.call(create_haplotypes, arg_list_)
        haps$n_haps()
    }, 5, info = paste("works with gene trees - ms-file"))

    expect_equal({
        arg_list_ <- c(list(haps_info = haps_gtrees(fn = .p("ms_out"))), arg_list)
        haps <- do.call(create_haplotypes, arg_list_)
        haps$n_haps()
    }, 5, info = "works with seg. sites - ms-file")

    expect_error(
        haps_gtrees(fn = .p("ms_out_err1")),
        regexp = "one or more chromosomes have no trees",
        info = "returns error when no gene trees provided - ms-file")

    expect_error(
        haps_ssites(fn = .p("ms_out_err1")),
        regexp = paste("One or more seg. sites matrices from a ms-style",
                       "file output have no haplotype information specified."),
        info = "returns error when no seg. sites provided - ms-file")


    expect_error({
        arg_list_ <- c(list(haps_info = haps_gtrees(fn = .p("ms_out_err3"))), arg_list)
        do.call(create_haplotypes, arg_list_)
    },
    regexp = paste("A coalescent string appears to include",
                   "recombination but does not include sizes for each region."),
    info = "returns error when improper gene trees provided - ms-file")


    expect_error({
        arg_list_ <- c(list(haps_info = haps_gtrees(fn = .p("ms_out_err3"))), arg_list)
        do.call(create_haplotypes, arg_list_)
    },
    regexp = paste("A coalescent string appears to include",
                   "recombination but does not include sizes for each region."),
    info = "returns error when improper gene trees provided - ms-file")


    expect_error({
        arg_list_ <- c(list(haps_info = haps_gtrees(fn = .p("ms_out_err2"))), arg_list)
        do.call(create_haplotypes, arg_list_)
    },
    regexp = "One or more trees have differing tip labels.",
    info = "returns error when improper gene trees provided - ms-file")


    expect_error(
        haps_ssites(fn = test_path(paste0("files/", "ms_out_err2.txt"))),
        regexp = paste("the listed positions for each site \\(line starting",
                       "with 'positions:'\\) does not have a length that's",
                       "the same as the \\# sites as given by the line",
                       "starting with 'segsites:'"),
        info = paste("returns error with improper seg. sites positions - ms-file"))


    expect_error(
        haps_ssites(fn = test_path("files/ms_out_err3.txt")),
        regexp = paste("the listed number of sites \\(line starting with",
                       "'segsites:'\\) does not agree with the number of",
                       "items in the 3th line of segregating sites",
                       "info \\(ones filled with 0s and 1s\\)"),
        info = paste("returns error with improper seg. sites matrices - ms file"))


    reference2 <- create_genome(3, 1000)

    expect_error({
        arg_list_ <- c(list(haps_info = haps_gtrees(fn = .p("ms_out"))), arg_list)
        arg_list_$reference <- reference2
        do.call(create_haplotypes, arg_list_)
    },
    regexp = paste("A coalescent string appears to include recombination",
                   "but the combined sizes of all regions don't match the",
                   "size of the chromosome"))

})



test_that("haplotype creation returns error with improper ref_genome input", {
    .p <- function(x) test_path(sprintf("files/%s.txt", x))
    expect_error({
        arg_list_ <- c(list(haps_info = haps_gtrees(fn = .p("ms_out"))), arg_list)
        arg_list_$reference <- list(1, 2)
        do.call(create_haplotypes, arg_list_)
    },
    regexp = paste("For the `create_haplotypes` function in jackalope,",
                   "argument `reference` must be a \"ref_genome\" object."))
})

test_that("seg. sites produces the correct number of mutations", {
    ms_file <- test_path("files/ms_out.txt")
    # Have to make bigger ref. genome to make sure site positions don't overlap.
    reference2 <- create_genome(3, 100e3)
    arg_list_ <- c(list(haps_info = haps_ssites(fn = ms_file)), arg_list)
    arg_list_$reference <- reference2
    haps <- do.call(create_haplotypes, arg_list_)

    msf <- readLines(ms_file)[-1:-2]
    msf <- msf[grepl("^0|^1", msf)]
    msf <- do.call(c, strsplit(msf, ""))
    n_muts_by_hap <-
        sapply(0:4, function(i) nrow(jackalope:::view_mutations(haps$ptr(), i)))
    expect_equal(sum(as.integer(msf)), sum(n_muts_by_hap))
})





# ==============================================================================`
# ==============================================================================`

# GENERAL NONSENSE -----

# ==============================================================================`
# ==============================================================================`


test_that("errors occur when nonsense is input to haps_gtrees", {

    expect_error(haps_gtrees(NULL, NULL),
                 regexp = "either argument `obj` or `fn` must be provided")
    expect_error(haps_gtrees(1, 1),
                 regexp = "only one argument \\(`obj` or `fn`\\) should be provided")
    expect_error(haps_gtrees(1, NULL),
                 regexp = paste("argument `obj` must be \\(1\\) a list with a",
                                "`trees` field present or \\(2\\) a list of lists,",
                                "each sub-list containing a `trees` field of length 1"))
    expect_error(haps_gtrees(NULL, 1),
                 regexp = "argument `fn` must be NULL or a single string")

})
test_that("errors occur when nonsense is input to haps_ssites", {

    expect_error(haps_ssites(NULL, NULL),
                 regexp = "either argument `obj` or `fn` must be provided")
    expect_error(haps_ssites(1, 1),
                 regexp = "only one argument \\(`obj` or `fn`\\) should be provided")
    expect_error(haps_ssites(1, NULL),
                 regexp = paste("argument `obj` must be NULL or a list with a",
                                "`seg_sites` field present"))
    expect_error(haps_ssites(NULL, 1),
                 regexp = "argument `fn` must be NULL or a single string")

})






# ==============================================================================`
# ==============================================================================`

# Writing gene trees -----

# ==============================================================================`
# ==============================================================================`


test_that("gene trees written properly by write_gtrees", {

    out_prefix <- paste0(tempdir(TRUE), "/trees")

    # Write gene trees to file based on known ms-style file:
    write_gtrees(haps_gtrees(fn = test_path("files/ms_out.txt")), out_prefix)

    # Newly written gene tree strings:
    wr_str <- strsplit(paste(readLines(paste0(out_prefix, ".trees")),
                             collapse = ""), "//")[[1]]
    wr_str <- wr_str[grepl("^\\[", wr_str)]

    # original file, split by chromosome:
    og_str <- strsplit(strsplit(paste(readLines(test_path("files/ms_out.txt"))[-1:-2],
                                      collapse = "\n"), "//\n")[[1]], "\n")[-1]
    # Now get just the gene trees:
    og_str <- sapply(og_str, function(x) paste(x[grepl("^\\[", x)], collapse = ""))

    expect_identical(wr_str, og_str)

})

