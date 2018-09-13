
# doc start ----
#' Create variants from a reference genome.
#'
#' @param reference A \code{ref_genome} object from which to generate variants.
#'     This argument is required.
#' @param method Method to use for generating variants.
#'     Options are as follows:
#'     \describe{
#'         \item{`"phylo"`}{a single phylogenetic tree from a `phylo` object.}
#'         \item{`"coal_obj"`}{a coalescent-simulator object from the `scrm` or `coala`
#'             packages.}
#'         \item{`"ms_file"`}{a file containing output from a coalescent simulator in the
#'             format of the `ms` program.}
#'         \item{`"newick"`}{a NEWICK file containing a phylogenetic tree.}
#'         \item{`"theta"`}{an estimate for theta, the population-scaled mutation rate.}
#'         \item{`"vcf"`}{a variant call format (VCF) file that directly specifies
#'             variants.
#'             \strong{\emph{NOTE:}} If this method is chosen, all arguments other than
#'             `reference`, `method`, and `method_args` are ignored.}
#'     }
#' @param method_args List of arguments used for the given method.
#'     See Details for which arguments are used for each method.
#' @param mevo_obj A `mevo` object that stores molecular-evolution information.
#'     See \code{\link{make_mevo}} for more information.
#' @param n_cores Number of cores to use for parallel processing.
#'     This argument is ignored if OpenMP is not enabled.
#'     Cores are spread across sequences, so it
#'     doesn't make sense to supply more cores than sequences in the reference genome.
#'     Defaults to `1`.
#' @param show_progress Boolean for whether to show a progress bar during processing.
#'     Defaults to `FALSE`.
#'
#'
#'
# doc end ----
create_variants <- function(reference,
                            method,
                            method_args,
                            mevo_obj,
                            n_cores = 1,
                            show_progress = FALSE) {

    method <- match.arg(method, c("phylo", "coal_obj", "ms_file", "newick", "theta",
                                  "vcf"))

    # ---------*
    # --- check types ----
    # ---------*

    if (!inherits(reference, "ref_genome")) {
        stop("\nCreating variants can only be done to a ref_genome object.",
             call. = FALSE)
    }
    ref_genome_ptr <- ref_genome$genome
    if (!inherits(ref_genome_ptr, "externalptr")) {
        stop("\nYou're attempting to create variants using a ref_genome object with ",
             "a genome field that is not an externalptr. ",
             "Restart by reading a FASTA file or by simulating a genome. ",
             "And do NOT change the genome field manually.",
             call. = FALSE)
    }
    if (!single_whole_number(n_cores, .min = 1)) {
        stop("\nThe n_cores argument supplied to create_variants is not a single",
             "whole number greater than or equal to 1.",
             call. = FALSE)
    }


    # ---------*
    # --- phylo methods ----
    # ---------*

    if (method != "vcf") {

        chunk_size <- mevo_obj$chunk_size

        phylo_info_ptr <- do.call(make_phylo_info,
                                  c(method_args, list(method = method,
                                                      chunk_size = chunk_size)))

        # -------+
        # Make sampler_base_ptr
        # -------+
        sampler_base_ptr <- mevo_obj$to_ptr()

        # -------+
        # Make Gamma matrices (for mutation-rate variability among sites):
        # -------+
        gamma_mats <- mevo_obj$gamma_mats

        # -------+
        # Make variants pointer:
        # -------+
        if (chunk_size > 0) {
            variant_ptr <- evolve_seqs_chunk(
                ref_genome_ptr,
                sampler_base_ptr,
                phylo_info_ptr,
                gamma_mats,
                n_cores,
                show_progress)
        } else {
            variant_ptr <- evolve_seqs(
                ref_genome_ptr,
                sampler_base_ptr,
                phylo_info_ptr,
                gamma_mats,
                n_cores,
                show_progress)
        }

    # ---------*
    # --- vcf method ----
    # ---------*

    } else {

        stop("\nVCF files not yet implemented in `create_variants`.", call. = TRUE)
        if (!inherits(vcf_file, "character") | length(vcf_file) != 1) {
            stop("\nThe vcf_file argument to create_variants must be a single string.",
                 call. = FALSE)
        }

    }

    var_obj <- variants$new(variants_ptr, ref_genome_ptr)

    return(var_obj)

}
# function end ----






