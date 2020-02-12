
#' Check PacBio arguments for validity.
#'
#'
#'
#' @noRd
#'
check_pacbio_args <- function(obj,
                              n_reads,
                              haplotype_probs,
                              sep_files,
                              compress,
                              comp_method,
                              n_threads,
                              read_pool_size,
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
                              custom_read_lengths,
                              prob_dup,
                              show_progress) {

    if (!inherits(obj, c("ref_genome", "haplotypes"))) {
        stop("\nWhen providing info for the PacBio sequencer, ",
             "the object providing the sequence information should be ",
             "of class \"ref_genome\" or \"haplotypes\".", call. = FALSE)
    }

    for (x in c("n_reads", "n_threads", "read_pool_size",
                "max_passes", "min_read_length")) {
        z <- eval(parse(text = x))
        if (!single_integer(z, 1)) err_msg("pacbio", x, "a single integer >= 1")
    }
    for (x in c("prob_thresh", "ins_prob", "del_prob", "sub_prob", "prob_dup")) {
        z <- eval(parse(text = x))
        if (!single_number(z, 0, 1)) err_msg("pacbio", x, "a single number in range [0,1].")
    }
    if (sum(c(ins_prob, del_prob, sub_prob)) > 1) {
        stop("\nWhen providing info for the PacBio sequencer, the insertion, ",
             "deletion, and substitution probabilities cannot sum to > 1.",
             call. = FALSE)
    }

    if (!is_type(chi2_params_s, "numeric", 5)) {
        err_msg("pacbio", "chi2_params_s", "a numeric vector of length 5.")
    }
    for (x in c("chi2_params_n", "lognorm_read_length")) {
        z <- eval(parse(text = x))
        if (!is_type(z, "numeric", 3)) err_msg("pacbio", x, "a numeric vector of length 3.")
    }
    for (x in c("sqrt_params", "norm_params")) {
        z <- eval(parse(text = x))
        if (!is_type(z, "numeric", 2)) err_msg("pacbio", x, "a numeric vector of length 2.")
    }

    if (!is.null(custom_read_lengths) && !is_type(custom_read_lengths, "matrix") &&
        !is_type(custom_read_lengths, "numeric")) {
        err_msg("pacbio", "custom_read_lengths", "NULL, a matrix, or a numeric vector.")
    }
    if (inherits(custom_read_lengths, "matrix") && ncol(custom_read_lengths) != 2) {
        stop("\nWhen providing info for the PacBio sequencer, ",
             "if the `custom_read_lengths` argument is a matrix, it should ",
             "have exactly 2 columns.", call. = FALSE)
    } else if (inherits(custom_read_lengths, "matrix")) {
        P <- custom_read_lengths[,2]
        if (any(P < 0) || all(P == 0)) {
            stop("\nWhen providing info for the PacBio sequencer, ",
                 "if the `custom_read_lengths` argument is a matrix, it should ",
                 "have exactly 2 columns, and the second column should contain ",
                 "no values < 0 and at least one value > 0.", call. = FALSE)
        }
    }

    if (!is.null(haplotype_probs) &&
        (!is_type(haplotype_probs, c("numeric", "integer")) ||
         any(haplotype_probs < 0) ||
         all(haplotype_probs == 0))) {
        err_msg("pacbio", "haplotype_probs", "NULL or a numeric/integer vector",
                "with no values < 0 and at least one value > 0")
    }
    if (!is_type(compress, "logical", 1) && !single_integer(compress, 1, 9)) {
        err_msg("pacbio", "compress", "a single logical or integer from 1 to 9")
    }
    if (!is_type(comp_method, "character", 1) || !comp_method %in% c("gzip", "bgzip")) {
        err_msg("pacbio", "comp_method", "\"gzip\" or \"bgzip\"")
    }
    if (!is_type(show_progress, "logical", 1)) {
        err_msg("pacbio", "show_progress", "a single logical")
    }
    if (!is_type(sep_files, "logical", 1)) {
        err_msg("pacbio", "sep_files", "a single logical")
    }
    # Checking proper haplotype_probs
    if (!is.null(haplotype_probs) && !inherits(obj, "haplotypes")) {
        stop("\nFor PacBio sequencing, it makes no sense to provide ",
             "a vector of probabilities of sequencing each haplotype if the ",
             "`obj` argument is of class \"ref_genome\". ",
             "Terminating here in case this was a mistake.", call. = FALSE)
    }
    if (!is.null(haplotype_probs) && inherits(obj, "haplotypes") &&
        length(haplotype_probs) != obj$n_haps()) {
        err_msg("pacbio", "haplotype_probs",
                "a vector of the same length as the number of haplotypes in the",
                "`obj` argument, if `obj` is of class \"haplotypes\".",
                "Use `obj$n_haps()` to see the number of haplotypes")
    }

    if (lognorm_read_length[1] < 0 || lognorm_read_length[3] < 0) {
        err_msg("pacbio", "lognorm_read_length",
                "a numeric vector of length 3 where items 1 and 3 cannot be < 0.")
    }

    invisible(NULL)
}







# pacbio doc -----
#' Create and write PacBio reads to FASTQ file(s).
#'
#'
#' From either a reference genome or set of variant haplotypes, create PacBio reads
#' and write them to FASTQ output file(s).
#' I encourage you to cite the reference below in addition to `jackalope` if you use
#' this function.
#'
#'
#' @section ID lines:
#' The ID lines for FASTQ files are formatted as such:
#'
#' `@<genome name>-<chromosome name>-<starting position>-<strand>`
#'
#' where `genome name` is always `REF` for reference genomes (as opposed to haplotypes).
#'
#'
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
#' @param sqrt_params  Vector containing the 2 parameters for the square root
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
#' @param prob_dup A single number indicating the probability of duplicates.
#'     Defaults to `0.0`.
#' @param read_pool_size The number of reads to store before writing to disk.
#'     Increasing this number should improve speed but take up more memory.
#'     Defaults to `100`.
#'
#' @return Nothing is returned.
#'
#' @export
#'
#'
#' @usage pacbio(obj,
#'        out_prefix,
#'        n_reads,
#'        chi2_params_s = c(0.01214, -5.12, 675, 48303.0732881,
#'                          1.4691051212330266),
#'        chi2_params_n = c(0.00189237136, 2.53944970, 5500),
#'        max_passes = 40,
#'        sqrt_params = c(0.5, 0.2247),
#'        norm_params = c(0, 0.2),
#'        prob_thresh = 0.2,
#'        ins_prob = 0.11,
#'        del_prob = 0.04,
#'        sub_prob = 0.01,
#'        min_read_length = 50,
#'        lognorm_read_length = c(0.200110276521, -10075.4363813,
#'                                17922.611306),
#'        custom_read_lengths = NULL,
#'        prob_dup = 0.0,
#'        haplotype_probs = NULL,
#'        sep_files = FALSE,
#'        compress = FALSE,
#'        comp_method = "bgzip",
#'        n_threads = 1L,
#'        read_pool_size = 100L,
#'        show_progress = FALSE,
#'        overwrite = FALSE)
#'
#'
#' @references
#' Stöcker, B. K., J. Köster, and S. Rahmann. 2016. SimLoRD: simulation of long
#' read data. \emph{Bioinformatics} \strong{32}:2704–2706.
#'
#' @examples
#' \donttest{
#' rg <- create_genome(10, 100e3, 100)
#' dir <- tempdir(TRUE)
#' pacbio(rg, paste0(dir, "/pacbio_reads"), n_reads = 100)
#' }
#'
pacbio <- function(obj,
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
                   prob_dup = 0.0,
                   haplotype_probs = NULL,
                   sep_files = FALSE,
                   compress = FALSE,
                   comp_method = "bgzip",
                   n_threads = 1L,
                   read_pool_size = 100L,
                   show_progress = FALSE,
                   overwrite = FALSE) {


    # Check for improper argument types:
    check_pacbio_args(obj, n_reads, haplotype_probs, sep_files,
                      compress, comp_method, n_threads, read_pool_size,
                      chi2_params_s, chi2_params_n, max_passes,
                      sqrt_params, norm_params,
                      prob_thresh, ins_prob, del_prob, sub_prob,
                      min_read_length, lognorm_read_length, custom_read_lengths,
                      prob_dup, show_progress)

    out_prefix <- path.expand(out_prefix)
    fns <- NULL
    # Doesn't make sense to have separate files for reference genome:
    if (inherits(obj, "ref_genome")) sep_files <- FALSE
    if (!sep_files) {
        fns <- paste0(out_prefix, "_R1.fq")
    } else {
        fns <- lapply(obj$hap_names(),
                      function(x) sprintf("%s_%s_R1.fq", out_prefix, x))
        fns <- c(fns, recursive = TRUE)
    }
    check_file_existence(fns, compress, overwrite)

    # Set compression level:
    if (is_type(compress, "logical", 1) && compress) compress <- 6 # default compression
    if (is_type(compress, "logical", 1) && !compress) compress <- 0 # no compression

    if (n_threads > 1 && compress > 0 && comp_method == "gzip") {
        stop("\nCompression using gzip cannot be performed using multiple threads. ",
             "Please use bgzip compression instead.")
    }

    if (!is.null(custom_read_lengths)) {
        if (inherits(custom_read_lengths, "matrix")) {
            read_lens <- custom_read_lengths[,1]
            read_probs <- custom_read_lengths[,2]
        } else {
            read_lens <- custom_read_lengths
            read_probs <- rep(1, length(custom_read_lengths))
        }
    } else {
        read_lens <- numeric(0)
        read_probs <- numeric(0)
    }

    if (is.null(haplotype_probs) && inherits(obj, "haplotypes")) {
        haplotype_probs <- rep(1, obj$n_haps())
    }

    # Assembling list of arguments for inner cpp function:
    args <- list(out_prefix = out_prefix,
                 sep_files = sep_files,
                 compress = compress,
                 comp_method = comp_method,
                 n_reads = n_reads,
                 n_threads = n_threads,
                 read_pool_size = read_pool_size,
                 chi2_params_s = chi2_params_s,
                 chi2_params_n = chi2_params_n,
                 max_passes = max_passes,
                 sqrt_params = sqrt_params,
                 norm_params = norm_params,
                 prob_thresh = prob_thresh,
                 prob_ins = ins_prob,
                 prob_del = del_prob,
                 prob_subst = sub_prob,
                 min_read_len = min_read_length,
                 scale = lognorm_read_length[3],
                 sigma = lognorm_read_length[1],
                 loc = lognorm_read_length[2],
                 read_lens = read_lens,
                 read_probs = read_probs,
                 prob_dup = prob_dup,
                 show_progress = show_progress)


    if (inherits(obj, "ref_genome")) {
        args <- c(args, list(ref_genome_ptr = obj$ptr()))
        args$sep_files <- NULL
        do.call(pacbio_ref_cpp, args)
    } else if (inherits(obj, "haplotypes")) {
        args <- c(args, list(hap_set_ptr = obj$ptr(),
                             haplotype_probs = haplotype_probs))
        do.call(pacbio_hap_cpp, args)
    } else {
        stop(paste("\nTrying to pass a `obj` argument to `pacbio` that's",
                   "not a \"ref_genome\" or \"haplotypes\" class."), call. = FALSE)
    }

    invisible(NULL)
}









