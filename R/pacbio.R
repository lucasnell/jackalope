
#' Check PacBio arguments for validity.
#'
#'
#'
#' @noRd
#'
check_pacbio_args <- function(seq_object,
                              n_reads,
                              variant_probs,
                              id_info,
                              compress,
                              n_cores,
                              read_chunk_size,
                              chi2_params_s,
                              chi2_params_n,
                              max_passes,
                              sqrt_params,
                              norm_params,
                              prob_thresh,
                              ins_prob,
                              del_prob,
                              sub_prob,
                              min_read_length,
                              lognorm_read_length,
                              custom_read_lengths) {

    err_msg <- function(par, ...) {
        stop(sprintf(paste("\nWhen providing info for the PacBio sequencer, the `%s`",
                           "argument must be %s."), par, paste(...)), call. = FALSE)
    }

    if (!inherits(seq_object, c("ref_genome", "variants"))) {
        stop("\nWhen providing info for the PacBio sequencer, ",
             "the object providing the sequence information should be ",
             "of class \"ref_genome\" or \"variants\".", call. = FALSE)
    }

    for (x in c("n_reads", "n_cores", "read_chunk_size",
                "max_passes", "min_read_length")) {
        z <- eval(parse(text = x))
        if (!single_integer(z, 1)) err_msg(x, "a single integer >= 1")
    }
    for (x in c("prob_thresh", "ins_prob", "del_prob", "sub_prob")) {
        z <- eval(parse(text = x))
        if (!single_number(z, 0, 1)) err_msg(x, "a single number in range [0,1].")
    }
    if (sum(c(ins_prob, del_prob, sub_prob)) > 1) {
        stop("\nWhen providing info for the PacBio sequencer, the insertion, ",
             "deletion, and substitution probabilities cannot sum to > 1.",
             call. = FALSE)
    }

    if (!is_type(chi2_params_s, "numeric", 5)) {
        err_msg("chi2_params_s", "a numeric vector of length 5.")
    }
    for (x in c("chi2_params_n", "lognorm_read_length")) {
        z <- eval(parse(text = x))
        if (!is_type(z, "numeric", 3)) err_msg(x, "a numeric vector of length 3.")
    }
    for (x in c("sqrt_params", "norm_params")) {
        z <- eval(parse(text = x))
        if (!is_type(z, "numeric", 2)) err_msg(x, "a numeric vector of length 2.")
    }

    if (!is.null(custom_read_lengths) && !is_type(custom_read_lengths, "matrix") &&
        !is_type(custom_read_lengths, "numeric")) {
        err_msg("custom_read_lengths", "NULL, a matrix, or a numeric vector.")
    }
    if (inherits(custom_read_lengths, "matrix") && ncol(custom_read_lengths) != 2) {
        stop("\nWhen providing info for the PacBio sequencer, ",
             "if the `custom_read_lengths` argument is a matrix, it should ",
             "have exactly 2 columns.", call. = FALSE)
    }

    if (!is.null(variant_probs) && !is_type(variant_probs, c("numeric", "integer"))) {
        err_msg("variant_probs", "NULL or a numeric/integer vector")
    }
    if (!is_type(id_info, "list")) err_msg("id_info", "a list.")
    if (!is_type(compress, "logical", 1)) err_msg("compress", "a single logical.")


    # Checking proper variant_probs
    if (!is.null(variant_probs) && !inherits(seq_object, "variants")) {
        stop("\nFor PacBio sequencing, it makes no sense to provide ",
             "a vector of probabilities of sequencing each variant if the ",
             "`seq_object` argument is of class \"ref_genome\". ",
             "Terminating here in case this was a mistake.", call. = FALSE)
    }
    if (!is.null(variant_probs) && inherits(seq_object, "variants") &&
        length(variant_probs) != seq_object$n_vars()) {
        err_msg("variant_probs",
                "a vector of the same length as the number of variants in the",
                "`seq_object` argument, if `seq_object` is of class \"variants\".",
                "Use `seq_object$n_vars()` to see the number of variants")
    }

    invisible(NULL)
}







# pacbio doc -----
#' Create and write PacBio reads to FASTQ file(s).
#'
#'
#' From either a reference genome or set of haploid variants, create PacBio reads
#' and write them to FASTQ output file(s).
#'
#' @inheritParams illumina
#' @param chi2_params_s  Vector containing the 5 parameters for the curve determining
#'     the scale parameter for the chi^2 distribution.
#'     Defaults to `c(0.01214, -5.12, 675, 48303.0732881, 1.4691051212330266)`.
#' @param chi2_params_n  Vector containing the 3 parameters for the function
#'     determining the n parameter for the chi^2 distribution.
#'     Defaults to `c(0.00189237136, 2.53944970, 5500)`.
#' @param max_passes  Maximal number of passes for one molecule.
#'     Defaults to `40`.
#' @param sqrt_params  Vector containing the 2 parameters for the sqare root
#'     function for the quality increase.
#'     Defaults to `c(0.5, 0.2247)`.
#' @param norm_params  Vector containing the 2 parameters for normal distributed
#'     noise added to quality increase square root function
#'     Defaults to `c(0, 0.2)`.
#' @param prob_thresh  Upper bound for the modified total error probability.
#'     Defaults to `0.2`.
#' @param ins_prob  Probability for insertions for reads with one pass.
#'     Defaults to `0.11`.
#' @param del_prob  Probability for deletions for reads with one pass.
#'     Defaults to `0.04`.
#' @param sub_prob  Probability for substitutions for reads with one pass.
#'     Defaults to `0.01`.
#' @param min_read_length  Minium read length for lognormal distribution.
#'     Defaults to `50`.
#' @param lognorm_read_length  Vector containing the 3 parameters for lognormal
#'     read length distribution.
#'     Defaults to `c(0.200110276521, -10075.4363813, 17922.611306)`.
#' @param custom_read_lengths  Sample read lengths from a vector or column in a
#'     matrix; if a matrix, the second column specifies the sampling weights.
#'     If `NULL`, it samples read lengths from the lognormal distribution
#'     using parameters in `lognorm_read_length`.
#'     Defaults to `NULL`.
#'
#' @return Nothing is returned.
#'
#' @export
#'
#' @examples
#'
#'
#'
pacbio <- function(seq_object,
                   out_prefix,
                   n_reads,
                   chi2_params_s = c(0.01214, -5.12, 675, 48303.0732881,
                                     1.4691051212330266),
                   chi2_params_n = c(0.00189237136, 2.53944970, 5500),
                   max_passes = 40,
                   sqrt_params = c(0.5, 0.2247),
                   norm_params = c(0, 0.2),
                   prob_thresh = 0.2,
                   ins_prob = 0.11,
                   del_prob = 0.04,
                   sub_prob = 0.01,
                   min_read_length = 50,
                   lognorm_read_length = c(0.200110276521, -10075.4363813, 17922.611306),
                   custom_read_lengths = NULL,
                   variant_probs = NULL,
                   id_info = list(),
                   compress = FALSE,
                   n_cores = 1L,
                   read_chunk_size = 100L) {

    out_prefix <- path.expand(out_prefix)

    # Check for improper argument types:
    check_pacbio_args(seq_object, n_reads, variant_probs, id_info,
                      compress, n_cores, read_chunk_size,
                      chi2_params_s, chi2_params_n, max_passes,
                      sqrt_params, norm_params,
                      prob_thresh, ins_prob, del_prob, sub_prob,
                      min_read_length, lognorm_read_length, custom_read_lengths)

    if (!is.null(custom_read_lengths)) {
        if (inherits(custom_read_lengths, "matrix")) {
            frag_lens <- custom_read_lengths[,1]
            frag_lens_probs <- custom_read_lengths[,2]
        } else {
            frag_lens <- custom_read_lengths
            frag_lens_probs <- rep(1, length(custom_read_lengths))
        }
    } else {
        frag_lens <- numeric(0)
        frag_lens_probs <- numeric(0)
    }

    if (is.null(variant_probs) && inherits(seq_object, "variants")) {
        variant_probs <- rep(1, seq_object$n_vars())
    }

    id_info <- do.call(id_line_info, id_info)


    # Assembling list of arguments for inner cpp function:
    args <- c(list(out_prefix = out_prefix,
                   compress = compress,
                   n_reads = n_reads,
                   n_cores = n_cores,
                   read_chunk_size = read_chunk_size,
                   chi2_params_s = chi2_params_s,
                   chi2_params_n = chi2_params_n,
                   max_passes = max_passes,
                   sqrt_params = sqrt_params,
                   norm_params = norm_params,
                   prob_thresh = prob_thresh,
                   ins_prob = ins_prob,
                   del_prob = del_prob,
                   sub_prob = sub_prob,
                   min_read_length = min_read_length,
                   lognorm_read_length = lognorm_read_length,
                   frag_lens = frag_lens,
                   frag_lens_probs = frag_lens_probs),
              id_info)

    if (inherits(seq_object, "ref_genome")) {
        args <- c(args, list(ref_genome_ptr = seq_object$genome))
        do.call(pacbio_ref_cpp, args)
    } else if (inherits(seq_object, "variants")) {
        args <- c(args, list(var_set_ptr = seq_object$genomes,
                             variant_probs = variant_probs))
        do.call(pacbio_var_cpp, args)
    } else {
        stop(paste("\nTrying to pass a `seq_object` argument to `pacbio` that's",
                   "not a \"ref_genome\" or \"variants\" class."), call. = FALSE)
    }

    invisible(NULL)
}









