
#' Convert to a XPtr<[Chunk]MutationSampler> object
#'
#' @noRd
#'
mevo_obj_to_ptr <- function(mevo_obj) {

    if (!single_integer(mevo_obj$chunk_size, 0)) {
        err_msg("create_variants", "mevo_obj", "a \"mevo\" object with a `chunk_size`",
                "field that's a single integer >= 0")
    }

    if (mevo_obj$chunk_size <= 0) {
        sampler_ptr <- make_mutation_sampler_base(mevo_obj$Q,
                                                  mevo_obj$pi_tcag,
                                                  mevo_obj$insertion_rates,
                                                  mevo_obj$deletion_rates)
    } else {
        sampler_ptr <- make_mutation_sampler_chunk_base(mevo_obj$Q,
                                                        mevo_obj$pi_tcag,
                                                        mevo_obj$insertion_rates,
                                                        mevo_obj$deletion_rates,
                                                        mevo_obj$chunk_size)
    }

    return(sampler_ptr)
}



# doc start ----
#' Create variants from a reference genome.
#'
#' @section Method arguments:
#' Below, I describe what the `method_info` should look like for each possible method.
#' \describe{
#'     \item{`method = "phylo"`}{One of the following object types is allowed:
#'         \itemize{
#'             \item A single \code{\link[ape]{phylo}} object that represents all
#'                 sequences in the genome.
#'             \item A `list` or `multiPhylo` object containing a `phylo` object for
#'                 each reference sequence.
#'                 Phylogenies will be assigned to sequences in the order provided.
#'         }
#'     }
#'     \item{`method = "coal_trees"`}{One of the following object types is allowed:
#'         \itemize{
#'             \item A single `list` with a `trees` field inside. This field must
#'                 contain a set of gene trees for each sequence.
#'             \item A list of lists, each sub-list containing a `trees` field of
#'                 length 1. The top-level list must be of the same length as the
#'                 number of sequences.
#'             \item A single string specifying the name of the file containing
#'                 the `ms`-style coalescent output with gene trees.
#'         }
#'         The top two options are designed after the `trees` fields in the output from
#'         the `scrm` and `coala` packages.
#'         (These packages are not required to be installed when installing
#'         `jackalope`.)
#'         To get gene trees, make sure to add `+ sumstat_trees()`
#'         to the `coalmodel` for `coala`, or
#'         make sure that `"-T"` is present in `args` for `scrm`.
#'         If using an output file from a command-line program like `ms`/`msms`,
#'         add the `-T` option.
#'     }
#'     \item{`method = "coal_sites"`}{One of the following object types is allowed:
#'         \itemize{
#'             \item A single `list` with a `seg_sites` field inside. This field must
#'                 contain a matrix for segregating sites for each sequence.
#'                 The matrix itself should contain the haplotype information, coded
#'                 using 0s and 1s: 0s indicate the ancestral state and 1s indicate
#'                 mutant.
#'                 The matrix column names should be numbers in the range (0,1) and
#'                 indicate the relative positions of the polymorphisms on the
#'                 chromosome.
#'             \item A single string specifying the name of the file containing
#'                 the `ms`-style coalescent output with segregating site info.
#'             \item A list containing `names` and `info` fields. The `info` field should
#'                 be one of the options above, and the `names` field provides a name
#'                 for each variant.
#'         }
#'         For what the `seg_sites` field should look like in a list, see output from the
#'         `scrm` or `coala` package.
#'         (These packages are not required to be installed when installing
#'         `jackalope`.)
#'     }
#'     \item{`method = "newick"`}{One or more string(s), each of which specifies
#'         a name of a NEWICK file containing a phylogeny.
#'         If one name is provided, that phylogeny will be used for all sequences.
#'         If more than one is provided, there must be a phylogeny for each sequence,
#'         and phylogenies will be assigned to sequences in the order provided.}
#'     \item{`method = "theta"`}{A named vector or list containing the fields `theta`
#'         and `n_vars`, specifying the theta parameter (population-scaled mutation rate)
#'         and number of desired variants, respectively.}
#'     \item{`method = "vcf"`}{Either (a) a single string specifying the name of
#'         the VCF file or (b) a list of arguments to pass to `vcfR::read.vcfR`.
#'         For the latter, the list can also contain the `print_chroms` field, which,
#'         if set to `TRUE`, prints all unique sequence names from the VCF file
#'         when VCF sequence names don't match those from the reference genome.
#'         This can be useful for troubleshooting.
#'         This method won't work if the package `vcfR` isn't installed.}
#' }
#'
#'
#' @param reference A \code{ref_genome} object from which to generate variants.
#'     This argument is required.
#' @param method Method to use for generating variants.
#'     Options are as follows:
#'     \describe{
#'         \item{`"phylo"`}{phylogenetic tree(s) from `phylo` object(s).}
#'         \item{`"coal_trees"`}{information from gene trees, either in the form of
#'             (1) coalescent-simulator object(s) from the `scrm` or `coala` package, or
#'             (2) a file containing output from a coalescent simulator in the
#'             format of the `ms` program.}
#'         \item{`"coal_sites"`}{information from matrices of segregating sites,
#'             either in the form of
#'             (1) coalescent-simulator object(s) from the `scrm` or `coala` package, or
#'             (2) a file containing output from a coalescent simulator in the
#'             format of the `ms` program.}
#'         \item{`"newick"`}{NEWICK file(s) containing phylogenetic tree(s).}
#'         \item{`"theta"`}{an estimate for theta, the population-scaled mutation rate.}
#'         \item{`"vcf"`}{a variant call format (VCF) file that directly specifies
#'             variants. This method does not work if the `vcfR` package isn't installed.
#'             \strong{\emph{NOTE:}} If this method is chosen, all arguments other than
#'             `reference`, `method`, and `method_info` are ignored.}
#'     }
#' @param method_info Object containing information used for the given method.
#'     See "Method arguments" section for which arguments are used for each method.
#' @param sub Output from one of the \code{\link{sub_models}} functions that organizes
#'     information for the substitution models.
#'     See `?sub_models` for more information on these models and
#'     their required parameters.
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
#' @param chunk_size The size of "chunks" of sequences to first sample uniformly
#'     before doing weighted sampling by rates for each sequence location.
#'     Uniformly sampling before doing weighted sampling dramatically speeds up
#'     the mutation process (especially for very long sequences) and has little
#'     effect on the sampling probabilities.
#'     Higher values will more closely resemble sampling without the uniform-sampling
#'     step, but will be slower.
#'     Set this to `0` to not uniformly sample first.
#'     From testing on a chromosome of length `1e6`, a `chunk_size` value of `100`
#'     offers a ~10x speed increase and doesn't differ significantly from sampling
#'     without the uniform-sampling step.
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
#' tree <- ape::rcoal(5)
#' v_phylo <- create_variants(r, "phylo", tree, sub_JC69(0.1))
#' v_theta <- create_variants(r, "theta", list(theta = 0.001, n_vars = 5),
#'                            sub_K80(0.1, 0.2))
#'
# doc end ----
create_variants <- function(reference,
                            method,
                            method_info,
                            sub = NULL,
                            ins = NULL,
                            del = NULL,
                            gamma_mats = NULL,
                            chunk_size = 100,
                            n_threads = 1,
                            show_progress = FALSE) {

    methods_ <- list(phylo = c("phylo", "coal_trees", "newick", "theta"),
                     non = c("coal_sites", "vcf"))

    method <- match.arg(method, as.character(do.call(c, methods_)))

    # ---------*
    # --- check types ----
    # ---------*

    if (!inherits(reference, "ref_genome")) {
        err_msg("create_variants", "reference", "a \"ref_genome\" object")
    }
    if (!inherits(reference$genome, "externalptr")) {
        err_msg("create_variants", "mevo_obj", "a \"mevo\" object with a `genome`",
                "field of class \"externalptr\".",
                "Restart by reading a FASTA file or by simulating a genome,",
                "and do NOT change the `genome` field manually")
    }

    # Do checks and organize molecular-evolution info into `mevo` object:
    mevo_obj <- create_mevo(reference, sub, ins, del, gamma_mats, chunk_size)

    if (!single_integer(n_threads, .min = 1)) {
        err_msg("create_variants", "n_threads", "a single integer >= 1")
    }
    if (!is_type(show_progress, "logical", 1)) {
        err_msg("create_variants", "show_progress", "a single logical")
    }
    # Check mevo_obj argument
    if (method != "vcf" && is.null(mevo_obj)) {
        err_msg("create_variants", "sub", "provided if you",
                "want to use a method other than \"vcf\".",
                "You should use one of the `sub_models` functions to create this object",
                "(see `?sub_models`)")
    }


    # ---------*
    # --- phylo methods ----
    # ---------*

    if (method %in% methods_$phylo) {

        # -------+
        # Make phylo_ptr
        # -------+
        phylo_info_ptr <- make_phylo_info(method, method_info,
                                          reference, mevo_obj)


        # -------+
        # Make sampler_base_ptr
        # -------+
        sampler_base_ptr <- mevo_obj_to_ptr(mevo_obj)

        # -------+
        # Make Gamma matrices (for mutation-rate variability among sites):
        # -------+
        gamma_mats <- mevo_obj$gamma_mats

        # -------+
        # Make variants pointer:
        # -------+
        if (mevo_obj$chunk_size > 0) {
            variants_ptr <- evolve_seqs_chunk(
                reference$genome,
                sampler_base_ptr,
                phylo_info_ptr,
                gamma_mats,
                n_threads,
                show_progress)
        } else {
            variants_ptr <- evolve_seqs(
                reference$genome,
                sampler_base_ptr,
                phylo_info_ptr,
                gamma_mats,
                n_threads,
                show_progress)
        }


    # ---------*
    # --- coal_sites method ----
    # ---------*
    } else if (method == "coal_sites") {

        variants_ptr <- read_coal_sites(method_info, reference, mevo_obj,
                                        n_threads, show_progress)

    # ---------*
    # --- vcf method ----
    # ---------*
    } else {

        variants_ptr <- read_vcf(reference, method_info)

    }

    var_obj <- variants$new(variants_ptr, reference$genome)

    return(var_obj)

}






