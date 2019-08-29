

# doc start ----
#' Create variants from a reference genome.
#'
#' Uses one of multiple methods to create haploid variants from a reference genome.
#' See \code{\link{vars_functions}} for the methods available.
#'
#'
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
#'     If you create a blank `variants` object, you can use its `add_vars` method
#'     to add variants manually.
#' @param sub Output from one of the \code{\link{sub_models}} functions that organizes
#'     information for the substitution models.
#'     See \code{\link{sub_models}} for more information on these models and
#'     their required parameters.
#'     This argument is ignored if you are using a VCF file to create variants.
#'     Passing `NULL` to this argument results in no substitutions.
#'     Defaults to `NULL`.
#' @param ins Output from the \code{\link{indels}} function that specifies rates
#'     of insertions by length.
#'     This argument is ignored if you are using a VCF file to create variants.
#'     Passing `NULL` to this argument results in no insertions.
#'     Defaults to `NULL`.
#' @param del Output from the \code{\link{indels}} function that specifies rates
#'     of deletions by length.
#'     This argument is ignored if you are using a VCF file to create variants.
#'     Passing `NULL` to this argument results in no deletions.
#'     Defaults to `NULL`.
#' @param epsilon Error control parameter for the "tau-leaping" approximation to
#'     the Doob–Gillespie algorithm, as used for the indel portion of the simulations.
#'     Smaller values result in a closer approximation.
#'     Larger values are less exact but faster.
#'     Values must be `> 0` and `< 1`.
#'     For more information, see the references below.
#'     Defaults to `0.03`.
#' @param n_threads Number of threads to use for parallel processing.
#'     This argument is ignored if OpenMP is not enabled.
#'     Threads are spread across chromosomes, so it
#'     doesn't make sense to supply more threads than chromosomes in the reference genome.
#'     Defaults to `1`.
#' @param show_progress Boolean for whether to show a progress bar during processing.
#'     Defaults to `FALSE`.
#'
#'
#' @export
#'
#'
#' @references
#' Doob, J. L. 1942. Topics in the theory of markoff chains.
#' \emph{Transactions of the American Mathematical Society} \strong{52}(1): 37–64.
#'
#' Gillespie, D. T. 1976. A general method for numerically simulating the stochastic time
#' evolution of coupled chemical reactions. \emph{Journal of Computational Physics}
#' \strong{22}(4): 403–434.
#'
#' Wieder, N., R. H. Fink, and F. von Wegner. 2011. Exact and approximate stochastic
#' simulation of intracellular calcium dynamics.
#' \emph{Journal of Biomedicine and Biotechnology} \strong{2011}: 572492.
#'
#' Cao, Y., D. T. Gillespie, and L. R. Petzold. 2006. Efficient step size
#' selection for the tau-leaping simulation method.
#' \emph{The Journal of Chemical Physics} \strong{124}(4): 044109.
#'
#'
#'
#' @examples
#' r <- create_genome(10, 1000)
#' v_phylo <- create_variants(r, vars_phylo(ape::rcoal(5)), sub_JC69(0.1))
#' v_theta <- create_variants(r, vars_theta(0.001, 5), sub_K80(0.1, 0.2))
#'
# doc end ----
create_variants <- function(reference,
                            vars_info,
                            sub = NULL,
                            ins = NULL,
                            del = NULL,
                            epsilon = 0.03,
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
        variants_ptr <- make_var_set(reference$ptr(), 0)
        var_obj <- variants$new(variants_ptr, reference$ptr())
        return(var_obj)
    }


    if (!inherits(reference$ptr(), "externalptr")) {
        err_msg("create_variants", "reference", "a \"ref_genome\" object with a `ptr`",
                "method that returns an object of class \"externalptr\".",
                "Restart by reading a FASTA file or by simulating a genome.")
    }
    if (!inherits(vars_info, do.call(c, vic))) {
        err_msg("create_variants", "vars_info", "NULL or one of the following classes:",
                paste(sprintf("\"%s\"", do.call(c, vic)), collapse = ", "))
    }
    # Check that sub, ins, or del info was passed if a non-VCF method is desired:
    vcf <- inherits(vars_info, vic$non[grepl("vcf", vic$non)])
    if (!vcf && is.null(sub) && is.null(ins) && is.null(del)) {
        stop("\nFor the `create_variants` function in jackalope, ",
             "if you are using a variant-creation method other than a VCF file, ",
             "you must provide input to the `sub`, `ins`, or `del` argument.",
             call. = FALSE)
    }

    # -----------------*
    # Do checks and organize molecular-evolution info
    # (or `NULL` if `sub` was not provided):
    # -----------------*
    if (!is.null(sub) && !is_type(sub, "sub_info")) {
        err_msg("create_variants", "sub", "NULL or a \"sub_info\" object")
    }
    if (!is.null(ins) && !is_type(ins, "indel_info")) {
        err_msg("create_variants", "ins", "NULL or a \"indel_info\" object")
    }
    if (!is.null(del) && !is_type(del, "indel_info")) {
        err_msg("create_variants", "del", "NULL or a \"indel_info\" object")
    }
    if (!single_number(epsilon) || epsilon <= 0 || epsilon >= 1) {
        err_msg("create_variants", "epsilon", "a single number > 0 and < 1")
    }
    # If sub is NULL and it's not a vcf file method, convert sub to rate-0 matrix:
    if (is.null(sub) && !vcf) sub <- sub_JC69(0)

    # Below will turn `NULL` into indel_info object with `numeric(0)` as `rates` method:
    if (is.null(ins)) ins <- indel_info$new(numeric(0))
    if (is.null(del)) del <- indel_info$new(numeric(0))


    if (!single_integer(n_threads, .min = 1)) {
        err_msg("create_variants", "n_threads", "a single integer >= 1")
    }
    if (!is_type(show_progress, "logical", 1)) {
        err_msg("create_variants", "show_progress", "a single logical")
    }

    # `to_var_set` is a method defined for each class of input for `vars_info`
    variants_ptr <- to_var_set(x = vars_info,
                               reference = reference,
                               sub = sub,
                               ins = ins,
                               del = del,
                               n_threads = n_threads,
                               show_progress = show_progress)

    var_obj <- variants$new(variants_ptr, reference$ptr())

    return(var_obj)

}




