
#' Create variants from a reference genome.
#'
#' @param reference A \code{ref_genome} object from which to generate variants.
#' @param n_variants The number of haploid variants to create.
#' @param phylo_control A list containing parameters for phylogenetic evolution.
#' @param vcf_file The file name for a VCF file to read from.
#' @param random_control A list containing parameters for random mutation addition.
#' @inheritParams create_genome
#'
#'
create_variants <- function(reference,
                            n_variants,
                            phylo_control,
                            vcf_file,
                            random_control,
                            n_cores = 1) {

    # ---------
    # Check for only one of phylo_control, vcf_file, and random_control being provided
    # ---------
    if ((!missing(phylo_control) & !missing(vcf_file)) |
        (!missing(phylo_control) & !missing(random_control)) |
        (!missing(random_control) & !missing(vcf_file))) {
        stop("\nIn create_variants, only one of the arguments phylo_control, vcf_file, ",
             "and random_control can be provided.",
             call. = FALSE)
    }
    if (missing(vcf_file) & missing(phylo_control) & missing(random_control)) {
        stop("\nIn create_variants, one of the following arguments must be provided: ",
             "phylo_control, vcf_file, or random_control.",
             call. = FALSE)
    }

    # ---------
    # Checking for proper types:
    # ---------
    if (!inherits(reference, "ref_genome")) {
        stop("\nCreating variants can only be done to a ref_genome object.",
             call. = FALSE)
    }
    reference_ptr <- ref_genome$genome
    if (!inherits(reference_ptr, "externalptr")) {
        stop("\nYou're attempting to create variants using a ref_genome object with ",
             "a genome field that is not an externalptr. ",
             "Restart by reading a FASTA file or by simulating a genome. ",
             "And do NOT change the genome field manually.",
             call. = FALSE)
    }
    if (!single_whole_number(n_variants, .min = 1)) {
        stop("\nThe n_variants argument supplied to create_variants is not a single",
             "whole number greater than or equal to 1.",
             call. = FALSE)
    }
    if (!single_whole_number(n_cores, .min = 1)) {
        stop("\nThe n_cores argument supplied to create_variants is not a single",
             "whole number greater than or equal to 1.",
             call. = FALSE)
    }


    # ---------
    # Running whatever method was chosen
    # ---------
    if (!missing(phylo_control)) {
        stop("\nPhylogenetic methods not yet implemented.", call. = FALSE)
    }
    if (!missing(vcf_file)) {
        stop("\nVCF file reading not yet implemented.", call. = FALSE)
    }
    if (!missing(random_control)) {
        variants_ptr <- random_variants(reference_ptr, n_variants, n_cores,
                                        random_control)
    }

    # ---------
    # Create output object
    # ---------

    var_obj <- variants$new(variants_ptr, reference_ptr)

    return(var_obj)

}





