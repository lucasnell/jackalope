# read profile fxns -----

#' Return table of information about built-in Illumina profiles.
#'
#' @noRd
#'
builtin_illumina_profiles <- function() {
    profiles <- read.csv(
        text = paste('"name","read_length","read","file_name","abbrev"',
                     '"Genome Analyzer I",36,1,"EmpR36R1","GA1"',
                     '"Genome Analyzer I",36,2,"EmpR36R2","GA1"',
                     '"Genome Analyzer I",44,1,"EmpR44R1","GA1"',
                     '"Genome Analyzer I",44,2,"EmpR44R2","GA1"',
                     '"Genome Analyzer II",50,1,"EmpR50R1","GA2"',
                     '"Genome Analyzer II",50,2,"EmpR50R2","GA2"',
                     '"MiniSeq TruSeq",50,1,"MiniSeqTruSeqL50","MinS"',
                     '"Genome Analyzer II",75,1,"EmpR75R1","GA2"',
                     '"Genome Analyzer II",75,2,"EmpR75R2","GA2"',
                     '"NextSeq 500 v2",75,1,"NextSeq500v2L75R1","NS50"',
                     '"NextSeq 500 v2",75,2,"NextSeq500v2L75R2","NS50"',
                     '"HiSeq 1000",100,1,"Emp100R1","HS10"',
                     '"HiSeq 1000",100,2,"Emp100R2","HS10"',
                     '"HiSeq 2000",100,1,"HiSeq2000L100R1","HS20"',
                     '"HiSeq 2000",100,2,"HiSeq2000L100R2","HS20"',
                     '"HiSeq 2500",125,1,"HiSeq2500L125R1","HS25"',
                     '"HiSeq 2500",125,2,"HiSeq2500L125R2","HS25"',
                     '"HiSeq 2500",150,1,"HiSeq2500L150R1filter","HS25"',
                     '"HiSeq 2500",150,2,"HiSeq2500L150R2filter","HS25"',
                     '"HiSeqX v2.5 PCR free",150,1,"HiSeqXPCRfreeL150R1","HSXn"',
                     '"HiSeqX v2.5 PCR free",150,2,"HiSeqXPCRfreeL150R2","HSXn"',
                     '"HiSeqX v2.5 TruSeq",150,1,"HiSeqXtruSeqL150R1","HSXt"',
                     '"HiSeqX v2.5 TruSeq",150,2,"HiSeqXtruSeqL150R2","HSXt"',
                     '"MiSeq v1",250,1,"EmpMiSeq250R1","MSv1"',
                     '"MiSeq v1",250,2,"EmpMiSeq250R2","MSv1"',
                     '"MiSeq v3",250,1,"MiSeqv3L250R1","MSv3"',
                     '"MiSeq v3",250,2,"MiSeqv3L250R2","MSv3"', sep = "\n"),
        stringsAsFactors = FALSE)
    profiles <- profiles[order(profiles$read_length),]
    return(profiles)
}




#' Find a sequencing system for a built-in profile file for a given read length.
#'
#' @inheritParams read_profile
#'
#' @noRd
#'
seq_sys_by_read_length <- function(read_length) {

    seq_sys <- ""
    if (read_length <= 44){
        seq_sys <- "GA1"
    } else if (read_length <= 75){
        seq_sys <- "GA2"
    } else if (read_length <= 100){
        seq_sys <- "HS20"
    } else if (read_length <= 150){
        seq_sys <- "HS25"
    } else if (read_length <= 250){
        seq_sys <- "MSv1"
    } else {
        stop("\nNo built-in Illumina profile can generate reads of length ",
             paste(read_length), ".", call. = FALSE)
    }

    return(seq_sys)
}




#' Find a built-in, ART-style Illumina profile file.
#'
#' @inheritParams read_profile
#'
#' @return The filename to be read from `read_profile`.
#'
#' @noRd
#'
find_profile_file <- function(seq_sys, read_length, read) {

    profile_df <- builtin_illumina_profiles()

    criteria1 <- profile_df$name == seq_sys | profile_df$abbrev == seq_sys
    if (sum(criteria1) == 0) {
        print(profile_df[!duplicated(profile_df$name), c("name", "abbrev")],
              row.names = FALSE)
        stop("\nThe desired Illumina platform name isn't available. ",
             "See printed data frame above for names and abbreviations of ",
             "those available.", call. = FALSE)
    }
    criteria2 <- criteria1 & profile_df$read_length >= read_length
    if (sum(criteria2) == 0) {
        cat(paste(unique(profile_df$read_length[criteria1]), collapse = " "), "\n")
        stop("\nFor the desired Illumina platform, this package doesn't have ",
             "a read length that's as long as you want. ",
             "See printed values above for lengths that are available.",
             call. = FALSE)
    }
    criteria3 <- criteria2 & profile_df$read == read
    if (sum(criteria3) == 0) {
        stop("\nFor the desired Illumina platform and read length, ",
             "this package only has a profile for read number ",
             paste(ifelse(read == 1, 2, 1)), ".", call. = FALSE)
    }

    profile_fn <- paste0(profile_df$file_name[criteria3][1], ".txt")
    profile_fn <- system.file("art_profiles", profile_fn, package = "gemino",
                              mustWork = TRUE)

    return(profile_fn)
}








#' Check profile's information for validity and format it for use in a sampler class.
#'
#' @param profile_info List of info read from the profile file.
#' @inheritParams read_profile
#'
#' @noRd
#'
format_profile <- function(profile_info, read_length) {

    pos <- sapply(profile_info, function(x) x$pos)
    if (min(pos) > 0) stop("\nMinimum profile position should be zero.", call. = FALSE)
    if (max(pos) < (read_length-1)) {
        stop("\nMaximum profile position should be >= read_length - 1.", call. = FALSE)
    }

    nts <- sapply(profile_info, function(x) x$nt)

    qual_probs <- rep(list(NA), 4)
    quals <- rep(list(NA), 4)
    names(qual_probs) <- names(quals) <- c("T", "C", "A", "G")

    for (nt in names(qual_probs)) {
        probs_nt <- lapply(profile_info[nts == nt], function(x) x$probs)
        qual_nt <- lapply(profile_info[nts == nt], function(x) x$quals)
        pos_nt <- sapply(profile_info[nts == nt], function(x) x$pos)
        if (!identical(0:(length(pos_nt) - 1), pos_nt)) {
            stop("\nFor nucleotide ", nt, " in the profile, the positions ",
                 "aren't a vector from 0 to length(positions) - 1.", call. = FALSE)
        }
        if (length(probs_nt) != length(qual_nt)) {
            stop("\nFor nucleotide ", nt, " in the profile, the number of ",
                 "positions for qualities isn't the same as for the quality ",
                 "probabilities.", call. = FALSE)
        }
        if (length(probs_nt) < read_length) {
            stop("\nFor nucleotide ", nt, " in the profile, it doesn't provide at ",
                 "least as many positions as your desired read length.", call. = FALSE)
        }
        if (any(sapply(probs_nt, length) != sapply(qual_nt, length))) {
            stop("\nFor nucleotide ", nt, " in the profile, at least ",
                 "one of the positions has a number of qualities that doesn't ",
                 "match with the number of quality probabilities.", call. = FALSE)
        }
        # Order by position:
        probs_nt <- probs_nt[order(pos_nt)]
        qual_nt <- qual_nt[order(pos_nt)]
        pos_nt <- pos_nt[order(pos_nt)]
        # Remove unnecessary positions if necessary:
        if (length(probs_nt) > read_length) {
            probs_nt <- probs_nt[1:read_length]
            qual_nt <- qual_nt[1:read_length]
            pos_nt <- pos_nt[1:read_length]
        }
        qual_probs[[nt]] <- probs_nt
        quals[[nt]] <- qual_nt
    }

    # Names no longer needed
    names(qual_probs) <- names(quals) <- NULL

    return(list(qual_probs = qual_probs, quals = quals))

}




#' Read an ART-style Illumina profile file.
#'
#' Note: this function will need to be run twice for paired-end reads.
#'
#' @param profile_fn File name of custom ART-style profile file for Illumina reads.
#' @param seq_sys The name of a built-in Illumina platform.
#' @param read_length Desired read length. Used to shorten the output vectors of
#'     qualities and quality-probabilities if `read_length` is less than the number of
#'     positions specified in the profile.
#' @param read Which read to get the profile from.
#'
#' @return A list of qualities and quality-probabilities for each nucleotide and
#'     position on read. I wouldn't recommend printing this list
#'     because it's quite lengthy.
#'
#' @noRd
#'
read_profile <- function(profile_fn, seq_sys, read_length, read) {

    if (!is.null(profile_fn) && !is.null(seq_sys)) {
        stop("\nFor Illumina sequencing, the user should never provide both a custom ",
             "profile file and a sequencing system.", call. = FALSE)
    }
    if (is.null(profile_fn) && is.null(seq_sys)) {
        seq_sys <- seq_sys_by_read_length(read_length)
    }
    if (is.null(profile_fn)) {
        profile_fn <- find_profile_file(seq_sys, read_length, read)
    }

    profile_str <- readLines(profile_fn)
    profile_str <- profile_str[grepl("^T|^C|^A|^G", profile_str)]
    profile_str <- strsplit(profile_str, "\t")
    profile_info <-
        lapply(seq(1, length(profile_str), 2),
               function(i) {
                   nn <- length(profile_str[[i]])
                   probs_ <- as.numeric(profile_str[[i+1]][3:nn])
                   if (nn > 3) probs_ <- probs_ - c(0, probs_[1:(nn-3)])
                   probs_ <- probs_ / sum(probs_)
                   quals_ <- as.integer(profile_str[[i]][3:nn])
                   if (!identical(profile_str[[i]][1:2], profile_str[[i+1]][1:2])) {
                       stop("\nInput profile file does not have proper format. ",
                            "The two lines specifying quality and distances should ",
                            "always have the same values for nucleotide and position.",
                            call. = FALSE)
                   }
                   if (length(profile_str[[i]]) != length(profile_str[[i+1]])) {
                       stop("\nInput profile file does not have proper format. ",
                            "The two lines specifying quality and distances should ",
                            "always have the same number of tab-delimited columns.",
                            call. = FALSE)
                   }
                   list(nt = profile_str[[i]][1],
                        pos = as.integer(profile_str[[i]][2]),
                        quals = quals_,
                        probs = probs_)
               })


    final_info <- format_profile(profile_info, read_length)

    return(final_info)

}




# info functions -----



#' Check Illumina arguments for validity.
#'
#'
#'
#' @noRd
#'
check_illumina_args <- function(seq_object, n_reads,
                                read_length, paired,
                                frag_mean, frag_sd,
                                matepair,
                                seq_sys, profile1, profile2,
                                ins_prob1, del_prob1,
                                ins_prob2, del_prob2,
                                frag_len_min, frag_len_max,
                                variant_probs, barcodes, prob_dup,
                                compress, n_cores, read_chunk_size,
                                show_progress) {

    # Checking types:

    if (!inherits(seq_object, c("ref_genome", "variants"))) {
        stop("\nWhen providing info for the Illumina sequencer, ",
             "the object providing the sequence information should be ",
             "of class \"ref_genome\" or \"variants\".", call. = FALSE)
    }

    for (x in c("read_length", "n_reads", "n_cores", "read_chunk_size")) {
        z <- eval(parse(text = x))
        if (!single_integer(z, 1)) err_msg("illumina", x, "a single integer >= 1")
    }
    for (x in c("paired", "matepair", "compress", "show_progress")) {
        z <- eval(parse(text = x))
        if (!is_type(z, "logical", 1)) err_msg("illumina", x, "a single logical.")
    }
    for (x in c("frag_mean", "frag_sd")) {
        z <- eval(parse(text = x))
        if (!single_number(z) || z <= 0) err_msg("illumina", x, "a single number > 0")
    }
    for (x in c("ins_prob1", "del_prob1", "ins_prob2", "del_prob2", "prob_dup")) {
        z <- eval(parse(text = x))
        if (!single_number(z, 0, 1)) err_msg("illumina", x, "a single number in range [0,1].")
    }
    for (x in c("seq_sys", "profile1", "profile2")) {
        z <- eval(parse(text = x))
        if (!is.null(z) && !is_type(z, "character", 1)) {
            err_msg("illumina", x, "NULL or a single string")
        }
    }
    for (x in c("frag_len_min", "frag_len_max")) {
        z <- eval(parse(text = x))
        if (!is.null(z) && !single_integer(z, 1)) {
            err_msg("illumina", x, "NULL or a single integer >= 1")
        }
    }
    if (!is.null(variant_probs) && !is_type(variant_probs, c("numeric", "integer"))) {
        err_msg("illumina", "variant_probs", "NULL or a numeric/integer vector")
    }
    if (!is.null(barcodes) && !is_type(barcodes, "character")) {
        err_msg("illumina", "barcodes", "NULL or a character vector")
    }

    # Checking for proper profile info:
    if (paired) {
        if (!is.null(profile1) && is.null(profile2)) {
            stop("\nFor Illumina paired-end reads, if you provide a custom profile for ",
                 "read 1, you must also provide a file for read 2.", call. = FALSE)
        } else if (is.null(profile1) && !is.null(profile2)) {
            stop("\nFor Illumina paired-end reads, if you provide a custom profile for ",
                 "read 2, you must also provide a file for read 1.", call. = FALSE)
        }
    } else if (!is.null(profile2)) {
        stop("\nFor Illumina single-end reads, it makes no sense to provide ",
             "a custom profile for read 2. ",
             "Terminating here in case this was a mistake.", call. = FALSE)
    }

    # Checking proper variant_probs
    if (!is.null(variant_probs) && !inherits(seq_object, "variants")) {
        stop("\nFor Illumina sequencing, it makes no sense to provide ",
             "a vector of probabilities of sequencing each variant if the ",
             "`seq_object` argument is of class \"ref_genome\". ",
             "Terminating here in case this was a mistake.", call. = FALSE)
    }
    if (!is.null(variant_probs) && inherits(seq_object, "variants") &&
        length(variant_probs) != seq_object$n_vars()) {
        err_msg("illumina", "variant_probs",
                "a vector of the same length as the number of variants in the",
                "`seq_object` argument, if `seq_object` is of class \"variants\".",
                "Use `seq_object$n_vars()` to see the number of variants")
    }
    # Similar checks for barcodes
    if (!is.null(barcodes)) {
        if (!inherits(seq_object, "variants") && length(barcodes) != 1) {
            stop("\nFor Illumina sequencing, it makes no sense to provide ",
                 "a vector of multiple barcodes if the `seq_object` argument is ",
                 "of class \"ref_genome\". ",
                 "Terminating here in case this was a mistake.", call. = FALSE)
        }
        if (inherits(seq_object, "variants") && length(barcodes) != seq_object$n_vars()) {
            err_msg("illumina", "barcodes",
                    "a vector of the same length as the number of variants in the",
                    "`seq_object` argument, if `seq_object` is of class \"variants\".",
                    "Use `seq_object$n_vars()` to see the number of variants")
        }
        n_weirdo_chars <- sapply(strsplit(barcodes, ""),
                                 function(x) sum(!x %in% c("T", "C", "A", "G")))
        if (any(n_weirdo_chars > 0)) {
            err_msg("illumina", "barcodes", "NULL or a character vector with only the characters",
                    "\"T\", \"C\", \"A\", and \"G\" present")
        }
    }

    invisible(NULL)
}








# Illumina doc -----
#' Create and write Illumina reads to FASTQ file(s).
#'
#'
#' From either a reference genome or set of haploid variants, create Illumina reads
#' from error profiles and write them to FASTQ output file(s).
#'
#' @section Sequencing profiles:
#' This section outlines how to use the `seq_sys`, `profile1`,
#' and `profile2` arguments.
#' If all arguments are `NULL` (their defaults), a sequencing system is chosen
#' based on the read length.
#' If, however, one or more arguments has been provided, then how they're provided
#' should depend on whether you want single- or paired-end reads.
#'
#' __For single-end reads__
#'
#' - `profile2` should be `NULL`.
#' - Only `seq_sys` or `profile1` should be provided, not both.
#'
#'  __For paired-end reads__
#'
#' - If providing `seq_sys`, don't provide either `profile1` or `profile2`.
#' - If providing `profile1`, you must also provide `profile2` (they can be the
#'   same if you want) and you cannot provide `seq_sys`.
#'
#'
#' @section Sequencing systems:
#' Sequencing system options are the following, where, for each system,
#' "name" is the full name, "abbrev" is the abbreviated name,
#' "max_len" indicates the maximum allowed read length,
#' and
#' "paired" indicates whether paired-end sequencing is allowed.
#'
#' | name                 | abbrev  | max_len | paired |
#' |:---------------------|:--------|:--------|:-------|
#' | Genome Analyzer I    | GA1     | 44      | Yes    |
#' | Genome Analyzer II   | GA2     | 75      | Yes    |
#' | HiSeq 1000           | HS10    | 100     | Yes    |
#' | HiSeq 2000           | HS20    | 100     | Yes    |
#' | HiSeq 2500           | HS25    | 150     | Yes    |
#' | HiSeqX v2.5 PCR free | HSXn    | 150     | Yes    |
#' | HiSeqX v2.5 TruSeq   | HSXt    | 150     | Yes    |
#' | MiniSeq TruSeq       | MinS    | 50      | No     |
#' | MiSeq v1             | MSv1    | 250     | Yes    |
#' | MiSeq v3             | MSv3    | 250     | Yes    |
#' | NextSeq 500 v2       | NS50    | 75      | Yes    |
#'
#'
#' @section ID line:
#' Possible values for the `id_line` argument are as follows:
#' \describe{
#' \item{instrument}{Instrument ID (string with the following characters allowed:
#'     a–z, A–Z, 0–9 and underscore). Defaults to `"SIM"`.}
#' \item{run_number}{Run number on instrument (numeric). Defaults to `1`.}
#' \item{flowcell_ID}{ID for flowcell (string with the following characters allowed:
#'     a–z, A–Z, 0–9). Defaults to `"FCX"`.}
#' \item{lane}{Lane number (numeric). Defaults to `1`.}
#' \item{tile}{Tile number (numeric). Defaults to `15`.}
#' \item{x_pos}{Run number on instrument (numeric). Defaults to `6329`.}
#' \item{y_pos}{X coordinate of cluster (numeric). Defaults to `1045`.}
#' \item{read}{Read number (numeric). Note that changing this doesn't do anything
#'     because it's set by the sequencing sampler as reads are created. Defaults to `1`.}
#' \item{is_filtered}{`"Y"` if the read is filtered (did not pass), `"N"` otherwise.
#'     Defaults to `"N"`.}
#' \item{control_number}{`0` when none of the control bits are on, otherwise it
#'     is an even number. On HiSeq X systems, control specification is not
#'     performed and this number is always `0`. Defaults to `0`.}
#' \item{sample_number}{Sample number from sample sheet (numeric). Defaults to `1`.}
#' }
#'
#'
#'
#' @param seq_object Sequencing object of class `ref_genome` or `variants`.
#' @param out_prefix Prefix for the output file(s), including entire path except
#'     for the file extension.
#' @param n_reads Number of reads you want to create.
#' @param read_length Length of reads.
#' @param paired Logical for whether to use paired-end reads.
#' @param frag_mean Mean of the Gamma distribution that generates fragment sizes.
#' @param frag_sd Standard deviation of the Gamma distribution that generates
#'     fragment sizes.
#' @param matepair Logical for whether to simulate mate-pair reads.
#'     This argument is ignored if `paired` is `FALSE`.
#'     Defaults to `FALSE`.
#' @param seq_sys Full or abbreviated name of sequencing system to use.
#'     See "Sequencing systems" section for options.
#'     See "Sequencing profiles" section for more information on how this argument,
#'     `profile1`, and `profile2` are used to specify profiles.
#'     Defaults to `NULL`.
#' @param profile1 Custom profile file for read 1.
#'     See "Sequencing profiles" section for more information on how this argument,
#'     `profile2`, and `seq_sys` are used to specify profiles.
#'     Defaults to `NULL`.
#' @param profile2 Custom profile file for read 2.
#'     See "Sequencing profiles" section for more information on how this argument,
#'     `profile1`, and `seq_sys` are used to specify profiles.
#'     Defaults to `NULL`.
#' @param ins_prob1 Insertion probability for read 1. Defaults to `0.00009`.
#' @param del_prob1 Deletion probability for read 1. Defaults to `0.00011`.
#' @param ins_prob2 Insertion probability for read 2. Defaults to `0.00015`.
#' @param del_prob2 Deletion probability for read 2. Defaults to `0.00023`.
#' @param frag_len_min Minimum fragment size. A `NULL` value results in the read length.
#'     Defaults to `NULL`.
#' @param frag_len_max Maximum fragment size.
#'     A `NULL` value results in `2^32-1`, the maximum allowed value.
#'     Defaults to `NULL`
#' @param variant_probs Relative probability of sampling each variant.
#'     This is ignored if sequencing a reference genome.
#'     `NULL` results in all having the same probability.
#'     Defaults to `NULL`.
#' @param barcodes Character vector of barcodes for each variant, or a single barcode
#'     if sequencing a reference genome. `NULL` results in no barcodes.
#'     Defaults to `NULL`.
#' @param prob_dup A single number indicating the probability of duplicates.
#'     Defaults to `0.02`.
#' @param compress A logical for whether to compress output FASTQ files using gzip.
#'     Defaults to `FALSE`.
#' @param n_cores The number of cores to use in processing.
#'     This argument is ignored if the package was not compiled with OpenMP.
#'     Defaults to `1`.
#' @param read_chunk_size The number of reads to store before writing to disk.
#'     Increasing this number should improve speed but take up more memory.
#'     Defaults to `100`.
#' @param show_progress Logical for whether to show a progress bar.
#'     Defaults to `FALSE`.
#'
#' @return Nothing is returned.
#'
#' @export
#'
#' @source
#' Huang, W., L. Li, J. R. Myers, and G. T. Marth. 2012. ART: a next-generation
#' sequencing read simulator. \emph{Bioinformatics} \strong{28}:593–594.
#'
#' @examples
#'
#'
#'
illumina <- function(seq_object,
                     out_prefix,
                     n_reads,
                     read_length,
                     paired,
                     frag_mean,
                     frag_sd,
                     matepair = FALSE,
                     seq_sys = NULL,
                     profile1 = NULL,
                     profile2 = NULL,
                     ins_prob1 = 0.00009,
                     del_prob1 = 0.00011,
                     ins_prob2 = 0.00015,
                     del_prob2 = 0.00023,
                     frag_len_min = NULL,
                     frag_len_max = NULL,
                     variant_probs = NULL,
                     barcodes = NULL,
                     prob_dup = 0.02,
                     compress = FALSE,
                     n_cores = 1L,
                     read_chunk_size = 100L,
                     show_progress = FALSE) {

    out_prefix <- path.expand(out_prefix)
    check_fastq(out_prefix, n_read_ends = ifelse(paired, 2, 1))

    # Check for improper argument types:
    check_illumina_args(seq_object, n_reads, read_length, paired,
                        frag_mean, frag_sd, matepair, seq_sys, profile1, profile2,
                        ins_prob1, del_prob1, ins_prob2, del_prob2,
                        frag_len_min, frag_len_max, variant_probs, barcodes, prob_dup,
                        compress, n_cores, read_chunk_size, show_progress)

    # Change mean and SD to shape and scale of Gamma distribution:
    frag_len_shape <- (frag_mean / frag_sd)^2
    frag_len_scale <- frag_sd^2 / frag_mean

    if (is.null(frag_len_min)) frag_len_min <- read_length
    # (Because it's coerced to an unsigned 32-bit integer, values >= 2^32 will
    #  make the integer "overflow" and change to something you don't want.)
    if (is.null(frag_len_max) || frag_len_max > 2^32 - 1) frag_len_max <- 2^32 - 1
    if (frag_len_min > frag_len_max) {
        stop("\nFragment length min can't be less than the max. ",
             "For computational reasons, both should also be < 2^32, ",
             "and if `frag_len_min` is not provided, it's automatically changed ",
             "to the read length.", call. = FALSE)
    }
    if (is.null(variant_probs) && inherits(seq_object, "variants")) {
        variant_probs <- rep(1, seq_object$n_vars())
    }
    if (is.null(barcodes) && inherits(seq_object, "variants")) {
        barcodes <- rep("", seq_object$n_vars())
    } else if (is.null(barcodes)) barcodes <- ""

    prof_info1 <- read_profile(profile1, seq_sys, read_length, 1)
    qual_probs1  <- prof_info1$qual_probs
    quals1 <- prof_info1$quals

    qual_probs2 <- NULL
    quals2 <- NULL
    if (paired) {
        prof_info2 <- read_profile(profile2, seq_sys, read_length, 2)
        qual_probs2 <- prof_info2$qual_probs
        quals2 <- prof_info2$quals
    } else {
        qual_probs2 <- list(list(numeric(0)))
        quals2 <- list(list(numeric(0)))
    }

    # Assembling list of arguments for inner cpp function:
    args <- list(out_prefix = out_prefix,
                 compress = compress,
                 n_reads = n_reads,
                 paired = paired,
                 matepair = matepair,
                 prob_dup = prob_dup,
                 n_cores = n_cores,
                 read_chunk_size = read_chunk_size,
                 frag_len_shape = frag_len_shape,
                 frag_len_scale = frag_len_scale,
                 frag_len_min = frag_len_min,
                 frag_len_max = frag_len_max,
                 qual_probs1 = as.name(quote(qual_probs1)),
                 quals1 = as.name(quote(quals1)),
                 ins_prob1 = ins_prob1,
                 del_prob1 = del_prob1,
                 qual_probs2 = as.name(quote(qual_probs2)),
                 quals2 = as.name(quote(quals2)),
                 ins_prob2 = ins_prob2,
                 del_prob2 = del_prob2,
                 barcodes = barcodes,
                 show_progress = show_progress)

    if (inherits(seq_object, "ref_genome")) {
        args <- c(args, list(ref_genome_ptr = seq_object$genome))
        do.call(illumina_ref_cpp, args)
    } else if (inherits(seq_object, "variants")) {
        args <- c(args, list(var_set_ptr = seq_object$genomes,
                             variant_probs = variant_probs))
        do.call(illumina_var_cpp, args)
    } else {
        err_msg("illumina", "`seq_object`", "a \"ref_genome\" or \"variants\" object")
    }

    invisible(NULL)
}









