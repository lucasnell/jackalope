
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


#' Checks for proper output and errors when creating variants from coalescent simulations
#'
#' coal_info can be a single string for a file or a coalescent object from scrm or coala.
#' If a file, you must manually change the file itself to get the proper errors below,
#' inputting different file names for errors vs proper output.
#'
#'
#'
test_coal_obj_sims <- function(coal_info, pkg, type) {

    # Tests this fxn can run, in order:
    tests <- c("run trees", "run sites", "error no trees", "error no sites",
               "error missing tree pos", "error missing tree tip", "error bad site pos")

    # Add path to filename if that's what's input
    if (inherits(coal_info, "character") && length(coal_info) == 1) {
        coal_info <- test_path(paste0("files/", coal_info))
    }

    if (type == tests[1]) {
        vars <- create_variants(reference, method = "coal_trees", coal_info, mevo_obj)
        expect_equal(vars$n_vars(), 5, info = paste("works with gene trees -", pkg))
    } else if (type == tests[2]) {
        vars <- create_variants(reference, method = "coal_sites", coal_info, mevo_obj)
        expect_equal(vars$n_vars(), 5, info = paste("works with seg. sites -", pkg))
    } else if (type == tests[3]) {
        if (inherits(coal_info, "list")) {
            coal_info$trees <- NULL
            err_msg_regexp <- paste("`method_info` must be \\(1\\) a list with a",
                                    "`trees` field present or \\(2\\) a list of",
                                    "lists, each sub-list containing a `trees`",
                                    "field of length 1.")
        } else err_msg_regexp <- "one or more sequences have no trees"
        expect_error(
            create_variants(reference, method = "coal_trees", coal_info, mevo_obj),
            regexp = err_msg_regexp,
            info = paste("returns error when no gene trees provided -", pkg))
    } else if (type == tests[4]) {
        if (inherits(coal_info, "list")) {
            coal_info$seg_sites <- NULL
            err_msg_regexp <- paste("argument `method_info` must be a list with",
                                    "a `seg_sites` field present")
        } else err_msg_regexp <- paste("One or more seg. sites matrices from a ms-style",
                                       "file output have no variant information specified.")
        expect_error(
            create_variants(reference, method = "coal_sites", coal_info, mevo_obj),
            regexp = err_msg_regexp,
            info = paste("returns error when no seg. sites provided -", pkg))
    } else if (type == tests[5]) {
        if (inherits(coal_info, "list")) {
            i <- which(sapply(coal_info$trees, length) > 1)[1]
            coal_info$trees[[i]][1] <- substr(coal_info$trees[[i]][1],
                                              regexpr("\\]", coal_info$trees[[i]][1])[[1]]+1,
                                              9999)
        }
        expect_error(
            create_variants(reference, method = "coal_trees", coal_info, mevo_obj),
            regexp = paste("A coalescent string appears to include",
                           "recombination but does not include sizes for each region."),
            info = paste("returns error when improper gene trees provided -", pkg))
    } else if (type == tests[6]) {
        if (inherits(coal_info, "list")) {
            i <- which(sapply(coal_info$trees, length) > 1)[1]
            strs <- coal_info$trees[[i]]
            pos <- as.integer(sapply(strs, function(x) substr(x, 2, regexpr("\\]", x)[[1]]-1)))
            tr <- ape::read.tree(text = strs)
            tr <- lapply(tr, ape::drop.tip, tip = "1")
            class(tr) <- "multiPhylo"
            coal_info$trees[[i]] <- sprintf("[%i]%s", pos, ape::write.tree(tr))
        }
        expect_error(
            create_variants(reference, method = "coal_trees", coal_info, mevo_obj),
            regexp = "the gene trees don't all have the same number of tips.",
            info = paste("returns error when improper gene trees provided -", pkg))
    } else if (type == tests[7]) {
        if (inherits(coal_info, "list")) {
            i <- which(sapply(coal_info$seg_sites, ncol) > 1)[1]
            mat <- coal_info$seg_sites[[i]]
            mat <- as.matrix(mat)
            colnames(mat) <- c("XX", colnames(mat)[-1])
            coal_info$seg_sites[[i]] <- mat
            expect_error(
                suppressWarnings({
                    create_variants(reference, method = "coal_sites", coal_info, mevo_obj)
                }),
                regexp = paste("Positions in one or more segregating-sites matrices",
                               "are producing NAs"),
                info = paste("returns error with improper seg. sites col names -", pkg))
        } else {
            expect_error(
                create_variants(reference, method = "coal_sites", coal_info, mevo_obj),
                regexp = paste("the listed positions for each site \\(line starting",
                               "with 'positions:'\\) does not have a length that's",
                               "the same as the \\# sites as given by the line",
                               "starting with 'segsites:'"),
                info = paste("returns error with improper seg. sites positions -", pkg))
        }
    } else stop("Unknown test type provided")

    invisible(NULL)
}










reference <- create_genome(3, 100)
mevo_obj <- create_mevo(reference, list(model = "JC69", lambda = 0.1),
                      list(rate = 0.1, max_length = 10),
                      list(rate = 0.1, max_length = 10))



test_that("variant creation works with ms-style file output", {

    pkg <- "ms-style output file"

    test_coal_obj_sims("ms_out.txt", pkg, "run trees")
    test_coal_obj_sims("ms_out.txt", pkg, "run sites")

    test_coal_obj_sims("ms_out_err1.txt", pkg, "error no trees")
    test_coal_obj_sims("ms_out_err1.txt", pkg, "error no sites")

    test_coal_obj_sims("ms_out_err3.txt", pkg, "error missing tree pos")
    test_coal_obj_sims("ms_out_err2.txt", pkg, "error missing tree tip")
    test_coal_obj_sims("ms_out_err2.txt", pkg, "error bad site pos")
    expect_error(
        create_variants(reference, method = "coal_sites",
                        test_path("files/ms_out_err3.txt"), mevo_obj),
        regexp = paste("the listed number of sites \\(line starting with",
                       "'segsites:'\\) does not agree with the number of",
                       "items in the 3th line of segregating sites",
                       "info \\(ones filled with 0s and 1s\\)"),
        info = paste("returns error with improper seg. sites matrices -", "ms file"))


    reference2 <- create_genome(3, 1000)
    mevo_obj2 <- create_mevo(reference2, list(model = "JC69", lambda = 0.1))

    expect_error(create_variants(reference2, method = "coal_trees",
                                 test_path("files/ms_out.txt"), mevo_obj2),
                 regexp = paste("A coalescent string appears to include recombination",
                                "but the combined sizes of all regions don't match the",
                                "size of the sequence"))

    test_that("variant creation returns error with improper ref_genome input", {
        expect_error(
            create_variants(list(1, 2), method = "coal_trees",
                            test_path("files/ms_out.txt"), mevo_obj),
            regexp = paste("For the `create_variants` function in jackalope,",
                           "argument `reference` must be a \"ref_genome\" object."))
    })

    test_that("seg. sites produces the correct number of mutations", {
        ms_file <- test_path("files/ms_out.txt")
        # Have to make bigger ref. genome to make sure site positions don't overlap.
        reference2 <- create_genome(3, 100e3)
        mevo_obj2 <- create_mevo(reference2, list(model = "JC69", lambda = 0.1))
        vars <- create_variants(reference2, method = "coal_sites", ms_file, mevo_obj)

        msf <- readLines(ms_file)[-1:-2]
        msf <- msf[grepl("^0|^1", msf)]
        msf <- do.call(c, strsplit(msf, ""))
        n_muts_by_var <-
            sapply(0:4, function(i) nrow(jackalope:::view_mutations(vars$genomes, i)))
        expect_equal(sum(as.integer(msf)), sum(n_muts_by_var))
    })

})





test_that("variant creation works with scrm coalescent object", {

    skip_if_not_installed("scrm")
    library(scrm)
    coal_obj <- scrm("5 3 -r 3.1 100 -t 10 -T -L")

    test_coal_obj_sims(coal_obj, "scrm", "run trees")
    test_coal_obj_sims(coal_obj, "scrm", "run sites")
    test_coal_obj_sims(coal_obj, "scrm", "error no trees")
    test_coal_obj_sims(coal_obj, "scrm", "error no sites")
    test_coal_obj_sims(coal_obj, "scrm", "error missing tree pos")
    test_coal_obj_sims(coal_obj, "scrm", "error missing tree tip")
    test_coal_obj_sims(coal_obj, "scrm", "error bad site pos")

})



test_that("variant creation works with coala coalescent object", {

    skip_if_not_installed("coala")
    library(coala)
    model <- coal_model(sample_size = 5, loci_number = 3, loci_length = 100) +
        feat_recombination(2) +
        feat_mutation(5) +
        sumstat_trees() +
        sumstat_seg_sites()
    coal_obj <- simulate(model)

    test_coal_obj_sims(coal_obj, "coala", "run trees")
    test_coal_obj_sims(coal_obj, "coala", "run sites")
    test_coal_obj_sims(coal_obj, "coala", "error no trees")
    test_coal_obj_sims(coal_obj, "coala", "error no sites")
    test_coal_obj_sims(coal_obj, "coala", "error missing tree pos")
    test_coal_obj_sims(coal_obj, "coala", "error missing tree tip")
    test_coal_obj_sims(coal_obj, "coala", "error bad site pos")

})

