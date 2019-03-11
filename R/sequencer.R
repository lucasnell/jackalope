
#' Use Mason defaults to generate default mismatch and quality info.
#'
#' @noRd
#'
illumina_mis_qual_defs <- function(read_length = 100) {

    # From Mason:
    prob_mismatch_scale <- 1.0
    prob_mismatch <- 0.004
    prob_mismatch_begin <- 0.002
    prob_mismatch_end <- 0.012
    position_raise <- 0.66
    mean_qual_begin <- 40
    mean_qual_end <- 39.5
    sd_qual_begin <- 0.05
    sd_qual_end <- 10
    mean_mismatch_qual_begin <- 39.5
    mean_mismatch_qual_end <- 30
    sd_mismatch_qual_begin <- 3
    sd_mismatch_qual_end <- 15

    # ---------*
    # Functions to calculate means and SDs by position:
    # ---------*
    by_pos_ <- function(x_begin, x_end, read_len) {
        x = 1:read_len / (read_len - 1)
        b = x_begin
        m = (x_end - x_begin)
        out <- (m * x + b)
        return(out)
    }
    # Piecewise version:
    by_pos_pw_ <- function(x_begin, x_avg, x_end, read_len, pos_raise) {

        # Compute probability at raise point.
        y_r  <- (2 * x_avg) - (pos_raise * x_begin) - x_end + (x_end * pos_raise)
        if (y_r < 0) {
            stop(paste("\nThe selected mismatch probabilities are not possible.",
                       "Trying to coerce the vector to have the desired average value",
                       "creates negative probabilites."),
                 call. = FALSE)
        }

        x = 1:read_len / (read_len - 1)
        b = ifelse(x < pos_raise, x_begin, y_r)
        m = ifelse(x < pos_raise, (y_r - x_begin) / pos_raise,
                   (x_end - y_r) / (1 - pos_raise))
        x = ifelse(x < pos_raise, x / pos_raise, x - pos_raise)

        out = m * x + b

        return(out)
    }

    # ---------*
    # Sequence-error info:
    # ---------*
    mis_probs <- by_pos_pw_(
        x_begin = prob_mismatch_begin, x_avg = prob_mismatch,
        x_end = prob_mismatch_end, read_len = read_length,
        pos_raise = position_raise)
    mis_probs <- mis_probs * prob_mismatch_scale


    # ---------*
    # Quality info
    # ---------*
    qual_means <- by_pos_(x_begin = mean_qual_begin, x_end = mean_qual_end,
                          read_len = read_length)
    qual_sds <- by_pos_(x_begin = sd_qual_begin, x_end = sd_qual_end,
                        read_len = read_length)
    qual_region_len <- 1

    mis_qual_means <- by_pos_(x_begin = mean_mismatch_qual_begin,
                              x_end = mean_mismatch_qual_end,
                              read_len = read_length)
    mis_qual_sds <- by_pos_(x_begin = sd_mismatch_qual_begin,
                            x_end = sd_mismatch_qual_end,
                            read_len = read_length)
    mis_qual_region_len <- 1


    output <- list(
        mis_probs = mis_probs,
        qual_means = qual_means,
        qual_sds = qual_sds,
        qual_region_len = qual_region_len,
        mis_qual_means = mis_qual_means,
        mis_qual_sds = mis_qual_sds,
        mis_qual_region_len = mis_qual_region_len)

    return(output)

}


#' Check sequencing info inputs for being proper types & wrap all but seq_type in a list.
#'
#' It also changes region lengths to 1 if passed as `NULL`.
#'
#'
#' @noRd
#'
check_seq_info_inputs <- function(seq_type,
                                  read_length,
                                  paired,
                                  frag_len_probs,
                                  frag_len_region_len,
                                  mis_probs,
                                  ins_probs,
                                  del_probs,
                                  error_region_len,
                                  ins_length_probs,
                                  del_length_probs,
                                  qual_means,
                                  qual_sds,
                                  qual_region_len,
                                  mis_qual_means,
                                  mis_qual_sds,
                                  mis_qual_region_len) {

    err_msg <- sprintf(paste("\nWhen providing info for the %s sequencer, the `%s`",
                             "argument must be %s."), seq_type, "%s", "%s")

    for (x in c("frag_len_probs",
                "mis_probs", "ins_probs", "del_probs",
                "ins_length_probs", "del_length_probs",
                "qual_means", "qual_sds",
                "mis_qual_means", "mis_qual_sds")) {
        z <- eval(parse(text = x))
        if (!is.null(z) && !is_type(z, "numeric")) {
            stop(sprintf(err_msg, x, "NULL or a numeric vector with no NAs"),
                 call. = FALSE)
        }
    }
    for (x in c("frag_len_region_len",
                "error_region_len",
                "qual_region_len",
                "mis_qual_region_len")) {
        z <- eval(parse(text = x))
        # If set to NULL, change to 1L
        if (is.null(z)) {
            assign(x, 1L)
        } else if (!single_integer(z, 1)) {
            stop(sprintf(err_msg, x, "NULL or a single whole number >= 1"), call. = FALSE)
        }
    }
    # Check that [mis_]qual_means and [mis_]qual_sds are the same size if both
    # are provided
    for (x in c("", "mis_")) {
        xs <- paste0(x, c("qual_means", "qual_sds"))
        zs <- lapply(xs, function(.x) eval(parse(text = .x)))
        if (!is.null(zs[[1]]) && !is.null(zs[[2]]) &&
            length(zs[[1]]) != length(zs[[2]])) {
            stop(sprintf(paste("\nWhen providing info for the %s sequencer, the",
                               "`%s` and `%s` arguments must be %s."),
                         seq_type, xs[1], xs[2], "the same length"),
                 call. = FALSE)
        }
    }

    if (seq_type == "Illumina") {

        if (!single_integer(read_length, 1)) {
            stop(sprintf(err_msg, "read_length", "a single whole number >= 1"), call. = FALSE)
        }
        if (!is_type(paired, "logical", 1)) {
            stop(sprintf(err_msg, "paired", "a single logical"), call. = FALSE)
        }

        req_size <- ceiling(read_length / frag_len_region_len)
        if (!is.null(frag_len_probs) && length(frag_len_probs) < req_size) {
            stop(sprintf(err_msg, "frag_len_probs",
                         paste("NULL or a numeric vector of length",
                               ">= `ceiling(read_length / frag_len_region_len)`.",
                               "If you really want fragments to never be >= the ",
                               "read length, add zeros to the end of `frag_len_probs`",
                               "until it's long enough.",
                               "(Note that `frag_len_region_len = 1` by default.)")),
                 call. = FALSE)
        }

        req_size <- ceiling(read_length / error_region_len)
        for (x in c("mis_probs", "ins_probs", "del_probs")) {
            z <- eval(parse(text = x))
            if (!is.null(z) && length(z) < req_size) {
                stop(sprintf(err_msg, x,
                             paste("NULL or a numeric vector of length",
                                   ">= `ceiling(read_length / error_region_len)`")),
                     call. = FALSE)
            }
        }
    }

    out_list <- list(read_length = read_length,
                     paired = paired,
                     frag_len_probs = frag_len_probs,
                     frag_len_region_len = frag_len_region_len,
                     mis_probs = mis_probs,
                     ins_probs = ins_probs,
                     del_probs = del_probs,
                     error_region_len = error_region_len,
                     ins_length_probs = ins_length_probs,
                     del_length_probs = del_length_probs,
                     qual_means = qual_means,
                     qual_sds = qual_sds,
                     qual_region_len = qual_region_len,
                     mis_qual_means = mis_qual_means,
                     mis_qual_sds = mis_qual_sds,
                     mis_qual_region_len = mis_qual_region_len)
    return(out_list)

}

