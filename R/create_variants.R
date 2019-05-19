

# doc start ----
#' Create variants from a reference genome.
#'
#' Uses one of multiple methods to create haploid variants from a reference genome.
#' See \code{\link{vars_functions}} for the methods available.
#'
#'
#' @param reference A \code{ref_genome} object from which to generate variants.
#'     This argument is required.
#' @param vars_info Output from one of the \code{\link{vars_functions}}.
#'     These functions organize higher-level information for use here.
#'     See \code{\link{vars_functions}} for brief descriptions and links to each method.
#'     If this argument is `NULL`, all arguments other than `reference` are ignored,
#'     and an empty `variants` object with no variants is returned.
#'     This is designed for use when you'd like to add mutations manually.
#'     If you create a blank `variants` object, you can= use its `add_vars` method
#'     to add variants manually.
#' @param sub Output from one of the \code{\link{sub_models}} functions that organizes
#'     information for the substitution models.
#'     See \code{\link{sub_models}} for more information on these models and
#'     their required parameters.
#'     This argument is only allowed to be missing if you are using a VCF file to create
#'     variants.
#'     Defaults to `NULL`.
#' @param ins Output from the \code{\link{indels}} function that specifies rates
#'     of insertions by length.
#'     Passing `NULL` to this argument results in no insertions.
#'     Defaults to `NULL`.
#' @param del Output from the \code{\link{indels}} function that specifies rates
#'     of deletions by length.
#'     Passing `NULL` to this argument results in no deletions.
#'     Defaults to `NULL`.
#' @param gamma_mats Output from the \code{\link{site_var}} function that specifies
#'     variability in mutation rates among sites (for both substitutions and indels).
#'     Passing `NULL` to this argument results in no variability among sites.
#'     Defaults to `NULL`.
#' @param region_size Size of regions to break genome into for sampling mutation
#'     locations.
#'     This causes Gamma regions to be split into smaller sections (obviously the
#'     Gamma values themselves are not changed).
#'     Doing this splitting is useful because sampling within a region is more
#'     computationally costly than sampling among regions.
#'     Higher numbers will result in lower memory usage but slower speed.
#'     Defaults to `100`.
#' @param n_threads Number of threads to use for parallel processing.
#'     This argument is ignored if OpenMP is not enabled.
#'     Threads are spread across sequences, so it
#'     doesn't make sense to supply more threads than sequences in the reference genome.
#'     Defaults to `1`.
#' @param show_progress Boolean for whether to show a progress bar during processing.
#'     Defaults to `FALSE`.
#'
#'
#' @export
#'
#' @examples
#' r <- create_genome(10, 1000)
#' v_phylo <- create_variants(r, vars_phylo(ape::rcoal(5)), sub_JC69(0.1))
#' v_theta <- create_variants(r, vars_theta(0.001, 5), sub_K80(0.1, 0.2))
#'
# doc end ----
create_variants <- function(reference,
                            vars_info,
                            sub,
                            ins = NULL,
                            del = NULL,
                            gamma_mats = NULL,
                            region_size = 100,
                            n_threads = 1,
                            show_progress = FALSE) {

    # `vars_info` classes:
    vic <- list(phylo = c("phylo", "gtrees", "theta"),
                              non = c("ssites", "vcf"))
    vic <- lapply(vic, function(x) paste0("vars_", x, "_info"))

    # ---------*
    # --- check types ----
    # ---------*

    if (!inherits(reference, "ref_genome")) {
        err_msg("create_variants", "reference", "a \"ref_genome\" object")
    }

    # Make empty `variants` object, ignoring everything other than `reference` argument:
    if (is.null(vars_info)) {
        variants_ptr <- make_var_set(reference$genome, 0)
        var_obj <- variants$new(variants_ptr, reference$genome)
        return(var_obj)
    }


    if (!inherits(reference$genome, "externalptr")) {
        err_msg("create_variants", "mevo_obj", "a \"mevo\" object with a `genome`",
                "field of class \"externalptr\".",
                "Restart by reading a FASTA file or by simulating a genome,",
                "and do NOT change the `genome` field manually")
    }
    if (!inherits(vars_info, do.call(c, vic))) {
        err_msg("create_variants", "vars_info", "NULL or one of the following classes:",
                paste(sprintf("\"%s\"", do.call(c, vic)), collapse = ", "))
    }
    # If you're using a VCF, change `sub` to `NULL` bc it's not used:
    if (inherits(vars_info, vic$non[grepl("vcf", vic$non)])) sub <- NULL
    # Check that sub info was passed if a non-VCF method is desired:
    if (missing(sub)) {
        err_msg("create_variants", "sub", "provided if you",
                "want to create variants using any method other than a VCF file.",
                "You should use one of the `sub_models` functions to create the `sub`",
                "argument object (see `?sub_models`)")
    }

    # Do checks and organize molecular-evolution info into `mevo` object
    # (or `NULL` if `sub` was not provided):
    mevo_obj <- create_mevo(reference, sub, ins, del, gamma_mats, region_size)

    if (!single_integer(n_threads, .min = 1)) {
        err_msg("create_variants", "n_threads", "a single integer >= 1")
    }
    if (!is_type(show_progress, "logical", 1)) {
        err_msg("create_variants", "show_progress", "a single logical")
    }

    # `to_var_set` is a method defined for each class of input for `vars_info`
    variants_ptr <- to_var_set(x = vars_info, reference = reference,
                               mevo_obj = mevo_obj, n_threads = n_threads,
                               show_progress = show_progress)

    var_obj <- variants$new(variants_ptr, reference$genome)

    return(var_obj)

}




