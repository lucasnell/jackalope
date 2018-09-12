
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
#'             variants.}
#'     }
#' @param method_params List of arguments used for the given method.
#'     See Details for which arguments are used for each method.
#' @param sub_model Character indicating which substitution mutation model to use.
#'     Takes one of the following options, with the required parameters in parentheses:
#'     \describe{
#'         \item{`"TN93"`}{ `pi_tcag`, `alpha_1`, `alpha_2`, `beta` }
#'         \item{`"JC69"`}{ `lambda` }
#'         \item{`"K80"`}{ `alpha`, `beta` }
#'         \item{`"F81"`}{ `pi_tcag` }
#'         \item{`"HKY85"`}{ `pi_tcag`, `alpha`, `beta` }
#'         \item{`"F84"`}{ `pi_tcag`, `beta`, `kappa` }
#'         \item{`"GTR"`}{ `pi_tcag`, `abcdef` }
#'         \item{`"UNREST"`}{ `Q` }
#'     }
#'
#'     Defaults to `"TN93"`.
#' @param sub_params A list containing the parameters for the specified
#'     substitution model.
#'     Defaults to `NULL`, which causes it to use all default parameters. See
#'     Details for default parameter values.
#' @param indel_params A list containing the parameters for indels.
#'     The following parameters are allowed:
#'     \describe{
#'         \item{`xi`}{Overall indel rate.}
#'         \item{`psi`}{Proportion of insertions to deletions.}
#'         \item{`rel_insertion_rates`}{Relative insertion rates.}
#'         \item{`rel_deletion_rates`}{Relative deletion rates.}
#'     }
#'     Defaults to `NULL`, which causes it to use all default parameters. See
#'     Details for default parameter values.
#' @param site_var_params List of parameters for generating variability in mutation
#'     rates among sites (for both substitutions and indels).
#'     A site's deviance from the average mutation rate is determined by its
#'     "gamma distance".
#'     A site's overall mutation rate is the mutation rate for that nucleotide
#'     (substitution + indel) multiplied by the site's gamma distance.
#'     There are two options for specifying gamma distances:
#'     \enumerate{
#'         \item Generate gamma distances from a Gamma distribution.
#'             This option requires the following arguments:
#'             \describe{
#'                 \item{`shape`}{Shape parameter for the Gamma distribution,
#'                     where the variance of the distribution is `1 / shape`.
#'                     The mean is fixed to 1.}
#'                 \item{`region_size`}{Size of regions where each site within that
#'                     region has the same gamma distance.}
#'             }
#'         \item Manually input matrices that specify the gamma distance and end points
#'             for regions each gamma distances refers to.
#'             This option requires the following argument:
#'             \describe{
#'                 \item{`mats`}{List of matrices, one for each sequence in the genome.
#'                     Each matrix should have two columns.
#'                     The first should contain the end points for each region.
#'                     The second should contain the gamma distances for each region.
#'                     Note that if gamma distances don't have a mean (weighted by
#'                     sequence length for each gamma distance) equal to 1, you're
#'                     essentially changing the overall mutation rate.}
#'             }
#'     }
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
#' @param n_cores Number of cores to use for parallel processing.
#'     This argument is ignored if OpenMP is not enabled.
#'     Cores are spread across sequences, so it
#'     doesn't make sense to supply more cores than sequences in the reference genome.
#'     Defaults to `1`.
#'
#' @param show_progress Boolean for whether to show a progress bar during processing.
#'     Defaults to `FALSE`.
#'
#' @inheritParams create_genome
#'
#'
# doc end ----
create_variants <- function(reference,
                            method,
                            method_params,
                            sub_model = "TN93",
                            sub_params = NULL,
                            indel_params = NULL,
                            site_var_params = NULL,
                            chunk_size = 100,
                            n_cores = 1,
                            show_progress = FALSE) {

    method <- match.arg(method, c("phylo", "coal_obj", "ms_file", "newick", "theta",
                                  "vcf"))
    sub_model <- match.arg(sub_model, c("TN93", "JC69", "K80", "F81", "HKY85", "F84",
                                        "GTR", "UNREST"))


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

        phylo_info_ptr <- do.call(make_phylo_info,
                                  c(method_args, list(method = method,
                                                      chunk_size = chunk_size)))

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Left off: Make sampler_base_ptr
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # sampler_base_ptr <- make_sampler_ptr(sub_params, indel_params, sub_model,
        #                                      chunk_size)

        # -------+
        # Make Gamma matrices (for mutation-rate variability among sites):
        # -------+
        seq_sizes <- view_ref_genome_seq_sizes(ref_genome_ptr)
        if (!is.null(site_var_params)) {
            if (all(c("shape", "region_size") %in% names(site_var_params))) {
                gamma_mats <- make_gamma_mats(seq_sizes,
                                              gamma_size_ = site_var_params$region_size,
                                              shape = site_var_params$shape)
            } else if ("mats" %in% names(site_var_params)) {

                err_msg <- paste("\nThe `mats` field inside the `site_var_params`",
                                "argument to the `create_variants` function needs to",
                                "be a list of matrices.")
                if (!inherits(site_var_params$mats, "list")) {
                    stop(err_msg, call. = FALSE)
                } else if (!all(sapply(site_var_params$mats, inherits,
                                       what = "matrix"))) {
                    stop(err_msg, call. = FALSE)
                }

                # Check matrices for proper end points and # columns:
                check_gamma_mats(site_var_params$mats, seq_sizes)

                gamma_mats <- site_var_params$mats

            } else {
                stop("\nThe `site_var_params` argument to `create_variants` ",
                     "must be a named list containing \"shape\" and \"region_size\" ",
                     "or \"mats\" as names.",
                     call. = FALSE)
            }
        } else {
            # This results in no variability among sites:
            gamma_mats <- make_gamma_mats(seq_sizes, gamma_size_ = 0, shape = 1)
        }

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






