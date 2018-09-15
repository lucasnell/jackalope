
# doc start ----
#' Create variants from a reference genome.
#'
#' @section Method arguments:
#' Below, I describe what the `method_info` should look like for each possible method.
#' \describe{
#'     \item{`method = "phylo"`}{A single \code{\link[ape]{phylo}} object.}
#'     \item{`method = "coal_obj"`}{A single `list` with a `trees` field inside.
#'         For what this `trees` field should look like, see output from the
#'         `scrm` or `coala` package. (These packages are not required to be installed
#'         when installing `gemino`.)}
#'     \item{`method = "ms_file"`}{A single string specifying the name of the file
#'         containing the `ms`-style coalescent output.}
#'     \item{`method = "newick"`}{A single string specifying the name of the NEWICK file
#'         containing the phylogeny.}
#'     \item{`method = "theta"`}{A named vector or list containing the fields `theta`
#'         and `n_vars`, specifying the theta parameter (population-scaled mutation rate)
#'         and number of desired variants, respectively.}
#'     \item{`method = "vcf"`}{A single string specifying the name of the VCF file.}
#' }
#'
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
#'             variants. This method does not work if the `vcfR` package isn't installed.
#'             \strong{\emph{NOTE:}} If this method is chosen, all arguments other than
#'             `reference`, `method`, and `method_info` are ignored.}
#'     }
#' @param method_info Object containing information used for the given method.
#'     See "Method arguments" section for which arguments are used for each method.
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
                            method_info,
                            mevo_obj,
                            n_cores = 1,
                            show_progress = FALSE) {

    methods_ <- list(phylo = c("phylo", "coal_obj", "ms_file", "newick", "theta"),
                     non = "vcf")

    method <- match.arg(method, do.call(c, methods_))

    # ---------*
    # --- check types ----
    # ---------*

    if (!inherits(reference, "ref_genome")) {
        stop("\nCreating variants can only be done to a ref_genome object.",
             call. = FALSE)
    }
    ref_genome_ptr <- ref_genome$genome
    if (!inherits(ref_genome_ptr, "externalptr")) {
        stop("\nYou're attempting to create variants using a \"ref_genome\" object with ",
             "a `genome` field that is not of class \"externalptr\". ",
             "Restart by reading a FASTA file or by simulating a genome. ",
             "And do NOT change the `genome` field manually.",
             call. = FALSE)
    }
    if (!single_whole_number(n_cores, .min = 1)) {
        stop("\nThe `n_cores` argument supplied to `create_variants` is not a single",
             "whole number greater than or equal to 1.",
             call. = FALSE)
    }


    # ---------*
    # --- phylo methods ----
    # ---------*

    if (method %in% methods_$phylo) {

        # -------+
        # Check mevo_obj argument
        # -------+
        mevo_obj_err <- FALSE
        if (missing(mevo_obj)) {
            mevo_obj_err <- TRUE
        } else if (!inherits(mevo_obj, "mevo")) {
            mevo_obj_err <- TRUE
        }
        if (mevo_obj_err) {
            stop("\nIf you want to use a method other than \"vcf\" in ",
                 "`create_variants`, you must provide the `mevo_obj` argument that is ",
                 "of class \"mevo\". ",
                 "You should use the `make_mevo` function to create this object.",
                 call. = FALSE)
        }

        # -------+
        # Make phylo_ptr
        # -------+
        seq_sizes <- reference$sizes
        n_seqs <- length(seq_sizes)
        phylo_info_ptr <- make_phylo_info(method, method_info,
                                          seq_sizes, n_seqs, mevo_obj$mu,
                                          mevo_obj$chunk_size)

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
        if (mevo_obj$chunk_size > 0) {
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
    # (It's the only non-phylogenetic method currently)
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






