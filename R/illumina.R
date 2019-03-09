
#'
#' You also need to organize higher-level stuff like...
#' - number of reads or coverage
#' - pooling, if any, to do
#' - % PCR duplicates (use `re_read` methods for this)
#' - barcodes for multiplexing and for which reads belong on which flow cell, lane, etc.
#'
#' See `SequenceIdentifierInfo` class inside `sequencer.h` to see the info required
#' for the FASTQ file identifier lines. (Much of this can probably be ignored.)
#'
#'
#' The classes `ReferenceIllumina` (typedef-ed from `Illumina_t<RefGenome>`)
#' and `VariantIllumina` are what you'll work with to do the sampling.
#' You want to focus on each class's `one_read` method.
#'
#' This is the most info you'd need for a sequencing object in `illumina.h`:
#'
#' For `VariantIllumina` only:
#'
#' - const VarSet& var_set
#' - const std::vector<double>& variant_probs
#'
#' For `ReferenceIllumina` only:
#'
#' - const RefGenome& seq_object
#'
#' For both, paired reads:
#'
#' - const double& frag_len_shape
#' - const double& frag_len_scale
#' - const uint32& frag_len_min_
#' - const uint32& frag_len_max_
#' - const std::vector<std::vector<std::vector<double>>>& qual_probs1
#' - const std::vector<std::vector<std::vector<uint8>>>& quals1
#' - const double& ins_prob1
#' - const double& del_prob1
#' - const std::vector<std::vector<std::vector<double>>>& qual_probs2
#' - const std::vector<std::vector<std::vector<uint8>>>& quals2
#' - const double& ins_prob2
#' - const double& del_prob2
#'
#' For both, single-end reads:
#'
#' - const double& frag_len_shape
#' - const double& frag_len_scale
#' - const uint32& frag_len_min_
#' - const uint32& frag_len_max_
#' - const std::vector<std::vector<std::vector<double>>>& qual_probs
#' - const std::vector<std::vector<std::vector<uint8>>>& quals
#' - const double& ins_prob
#' - const double& del_prob
#'
#'


#'
# make_illumina_info <- function(read_length = 100,
#                                paired = FALSE,
#                                frag_len_probs = NULL,
#                                frag_len_region_len = NULL,
#                                qual_probs = NULL,
#                                ins_probs = NULL,
#                                del_probs = NULL,
#                                error_region_len = NULL,
#                                ins_length_probs = NULL,
#                                del_length_probs = NULL,
#                                qual_means = NULL,
#                                qual_sds = NULL,
#                                qual_region_len = NULL,
#                                mis_qual_means = NULL,
#                                mis_qual_sds = NULL,
#                                mis_qual_region_len = NULL) {

source(textConnection(readLines("R/util.R")[8:30]))
source(textConnection(readLines("R/read_write.R")[240:452]))

#' read_length,
#' paired,
#' frag_mean,
#' frag_sd,
#' profile1,
#' profile2 = NULL,
#' ins_prob1 = 0.00009
#' del_prob1 = 0.00011
#' ins_prob2 = 0.00015
#' del_prob2 = 0.00023
#' frag_len_min = NULL,
#' frag_len_max = NULL

err_msg <- sprintf(paste("\nWhen providing info for the Illumina sequencer,",
                         "the `%s` argument must be %s."), "%s", "%s")

if (!single_integer(read_length, 1)) {
    stop(sprintf(err_msg, "read_length", "a single whole number >= 1"), call. = FALSE)
}
for (x in c("frag_len_min", "frag_len_max")) {
    z <- eval(parse(text = x))
    if (!is.null(z) && !single_integer(z, 1)) {
        stop(sprintf(err_msg, x, "NULL or a single whole number >= 1"), call. = FALSE)
    }
}
for (x in c("frag_mean",
            "frag_sd",
            "ins_prob1",
            "del_prob1",
            "ins_prob2",
            "del_prob2")) {
    z <- eval(parse(text = x))
    if (!is.null(z) && !single_number(z, 0)) {
        stop(sprintf(err_msg, x, "NULL or a single number >= 0"), call. = FALSE)
    }
}

if ((!is.null(frag_mean) && frag_mean <= 0) || (!is.null(frag_sd) && frag_sd <= 0)) {
    stop("\nFragment size mean and SD must both be > 0.", call. = FALSE)
}

frag_len_shape <- (frag_mean / frag_sd)^2
frag_len_scale <- frag_sd^2 / frag_mean

if (is.null(frag_len_min)) frag_len_min <- read_length
if (is.null(frag_len_max) || frag_len_max > 2^32 - 1) frag_len_max <- 2^32 - 1
if (frag_len_min > frag_len_max) {
    stop("\nFragment length min can't be less than the max.", call. = FALSE)
}
#' seq_system = NULL,
#' profile1 = NULL,
#' profile2 = NULL,

# LEFT OFF -----
if (!is_type(profile1, "character", 1)) {
    stop(sprintf(err_msg, "profile1", "a single string."), call. = FALSE)
}
if (!is.null(profile2) && !is_type(profile2, "character", 1)) {
    stop(sprintf(err_msg, "profile2", "NULL or a single string."), call. = FALSE)
}

if (paired && is.null(profile2)) {
    profile2 <- profile1
}

