
#' Use Mason defaults to generate sequencing information.
#'
#' Takes care of everything but the following fields:
#'
#' ```cpp
#' const T& seq_object
#' const std::vector<double>& frag_len_probs
#' const uint32& frag_len_region_len
#' const bool& paired
#' ```
#'
#' @noRd
#'
Mason_illumina_info <- function(read_length = 100,
                                prob_insert = 0.001,
                                prob_delete = 0.001,
                                prob_mismatch_scale = 1.0,
                                prob_mismatch = 0.004,
                                prob_mismatch_begin = 0.002,
                                prob_mismatch_end = 0.012,
                                position_raise = 0.66,
                                mean_qual_begin = 40,
                                mean_qual_end = 39.5,
                                sd_qual_begin = 0.05,
                                sd_qual_end = 10,
                                mean_mismatch_qual_begin = 39.5,
                                mean_mismatch_qual_end = 30,
                                sd_mismatch_qual_begin = 3,
                                sd_mismatch_qual_end = 15) {

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
    seq_error_info <- list()
    seq_error_info$mis_probs <- by_pos_pw_(
        x_begin = prob_mismatch_begin, x_avg = prob_mismatch,
        x_end = prob_mismatch_end, read_len = read_length,
        pos_raise = position_raise)
    seq_error_info$mis_probs <- seq_error_info$mis_probs * prob_mismatch_scale

    seq_error_info$ins_probs <- rep(prob_insert, read_length)
    seq_error_info$del_probs <- rep(prob_delete, read_length)
    seq_error_info$region_len <- 1

    seq_error_info$match_probs <- 1 - seq_error_info$mis_probs -
        seq_error_info$ins_probs - seq_error_info$del_probs

    # Indels are always of size 1
    ins_length_probs <- 1
    del_length_probs <- 1

    # ---------*
    # Quality info
    # ---------*
    qual_info <- list()
    qual_info$means <- by_pos_(x_begin = mean_qual_begin, x_end = mean_qual_end,
                               read_len = read_length)
    qual_info$sds <- by_pos_(x_begin = sd_qual_begin, x_end = sd_qual_end,
                             read_len = read_length)
    qual_info$region_len <- 1

    mis_qual_info <- list()
    mis_qual_info$means <- by_pos_(x_begin = mean_mismatch_qual_begin,
                                   x_end = mean_mismatch_qual_end,
                                   read_len = read_length)
    mis_qual_info$sds <- by_pos_(x_begin = sd_mismatch_qual_begin,
                                 x_end = sd_mismatch_qual_end,
                                 read_len = read_length)
    mis_qual_info$region_len <- 1


    output <- list(
        seq_error_info = seq_error_info,
        ins_length_probs = ins_length_probs,
        del_length_probs = del_length_probs,
        qual_info = qual_info,
        mis_qual_info = mis_qual_info,
        read_length = read_length)

    return(output)

}




# LEFT OFF -----------
#'
#' The next thing you need to figure out is how to input sequencing parameters to
#' the C++ code.
#' In `Mason_illumina_info` above, it creates much of the necessary info for Illumina
#' reads using default parameters in Mason.
#' You need to make a function that can use the Mason parameterization or can take
#' as inputs arguments that allow total control.
#' You need to do this for both Illumina and long-read sequencing (Nanopore and PacBio).
#'
#' You also need to provide a way for the user to control the fragmentation
#' process.
#' (Fragment sizes = read sizes in long-read sequencing.)
#'
#' You also need to organize higher-level stuff like barcodes for multiplexing and
#' for which reads belong on which flow cell, lane, etc.
#' See `SequenceIdentifierInfo` class inside `sequencer.h` to see the info required
#' for the FASTQ file identifier lines. (Much of this can probably be ignored.)
#'
#'
#' The template classes `IlluminaWGS_t` and `LongReadWGS_t` are what you'll work with
#' to do the sampling.
#' They've each been `typedef`ed for `RefGenome` and `VarGenome` classes, resulting in
#' `VariantIlluminaWGS` and `ReferenceIlluminaWGS` for `IlluminaWGS_t`, and
#' `VariantLongReadWGS` and `ReferenceLongReadWGS` for `LongReadWGS_t`.
#' You want to focus on each class's `one_read` method.
#'
#' This is the most info you'd need for a sequencing object in `sequencer.h`
#'   (last two not needed for long-reads):
#'
#'     const T& seq_object,
#'     const std::vector<double>& frag_len_probs,
#'     const uint32& frag_len_region_len,
#'     const List& seq_error_info,
#'         In this, there's...
#'         - const std::vector<double>& match_probs(seq_error_info["match_probs"]);
#'         - const std::vector<double>& mis_probs(seq_error_info["mis_probs"]);
#'         - const std::vector<double>& ins_probs(seq_error_info["ins_probs"]);
#'         - const std::vector<double>& del_probs(seq_error_info["del_probs"]);
#'         - const uint32& region_len_(seq_error_info["region_len"]);
#'     const std::vector<double>& ins_length_probs,
#'     const std::vector<double>& del_length_probs,
#'     const List& qual_info,
#'         In this, there's...
#'         - const std::vector<double>& means_(qual_info["means"]);
#'         - const std::vector<double>& sds_(qual_info["sds"]);
#'         - const uint32& region_len_(qual_info["region_len"]);
#'     const List& mis_qual_info,
#'         Same as above
#'     const uint32& read_length_,
#'     const bool& paired
#'

