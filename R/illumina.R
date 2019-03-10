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
                   probs_ <- probs_ - c(0, probs_[1:(nn-3)])
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
#' @inheritParams make_illumina_sampler
#'
#'
#' @noRd
#'
check_illumina_args <- function(seq_object, n_reads,
                                read_length, paired,
                                frag_mean, frag_sd,
                                seq_sys, profile1, profile2,
                                ins_prob1, del_prob1,
                                ins_prob2, del_prob2,
                                frag_len_min, frag_len_max,
                                variant_probs, barcodes, pcr_dups) {

    # Checking types:
    err_msg <- sprintf(paste("\nWhen providing info for the Illumina sequencer,",
                             "the `%s` argument must be %s."), "%s", "%s")

    if (!inherits(seq_object, c("ref_genome", "variants"))) {
        stop("\nWhen providing info for the Illumina sequencer, ",
             "the object providing the sequence information should be ",
             "of class \"ref_genome\" or \"variants\".", call. = FALSE)
    }

    for (x in c("read_length", "n_reads")) {
        z <- eval(parse(text = x))
        if (!single_integer(z, 1)) {
            stop(sprintf(err_msg, x, "a single integer >= 1"), call. = FALSE)
        }
    }
    if (!is_type(paired, "logical", 1)) {
        stop(sprintf(err_msg, "paired", "a single logical"), call. = FALSE)
    }
    for (x in c("frag_mean", "frag_sd", "ins_prob1", "del_prob1",
                "ins_prob2", "del_prob2")) {
        z <- eval(parse(text = x))
        if (!single_number(z) || z <= 0) {
            stop(sprintf(err_msg, x, "a single number > 0"), call. = FALSE)
        }
    }
    for (x in c("seq_sys", "profile1", "profile2")) {
        z <- eval(parse(text = x))
        if (!is.null(z) && !is_type(z, "character", 1)) {
            stop(sprintf(err_msg, x, "NULL or a single string"), call. = FALSE)
        }
    }
    for (x in c("frag_len_min", "frag_len_max")) {
        z <- eval(parse(text = x))
        if (!is.null(z) && !single_integer(z, 1)) {
            stop(sprintf(err_msg, x, "NULL or a single integer >= 1"), call. = FALSE)
        }
    }
    if (!is.null(variant_probs) && !is_type(variant_probs, c("numeric", "integer"))) {
        stop(sprintf(err_msg, "variant_probs", "NULL or a numeric/integer vector"),
             call. = FALSE)
    }
    if (!is.null(barcodes) && !is_type(barcodes, "character")) {
        stop(sprintf(err_msg, "barcodes", "NULL or a character vector"), call. = FALSE)
    }
    if (!single_number(pcr_dups, 0)) {
        stop(sprintf(err_msg, "pcr_dups", "a single number >= 0"), call. = FALSE)
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
        stop(sprintf(err_msg, "variant_probs",
                     paste("a vector of the same length as the number of variants",
                           "in the `seq_object` argument, if `seq_object` is",
                           "of class \"variants\". Use `seq_object$n_vars()`",
                           "to see the number of variants")), call. = FALSE)
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
            stop(sprintf(err_msg, "barcodes",
                         paste("a vector of the same length as the number of variants",
                               "in the `seq_object` argument, if `seq_object` is",
                               "of class \"variants\". Use `seq_object$n_vars()`",
                               "to see the number of variants")), call. = FALSE)
        }
        n_weirdo_chars <- sapply(strsplit(barcodes, ""),
                                 function(x) sum(!x %in% c("T", "C", "A", "G")))
        if (any(n_weirdo_chars > 0)) {
            stop(sprintf(err_msg, "barcodes",
                         paste("NULL or a character vector with only the characters",
                               "\"T\", \"C\", \"A\", and \"G\" present")),
                 call. = FALSE)
        }
    }

    invisible(NULL)
}



#' Organize arguments and return the pointer to the sampler object.
#'
#' @inheritParams sequence
#'
#' @noRd
#'
make_illumina_sampler <- function(seq_object, n_reads,
                                  read_length, paired,
                                  frag_mean, frag_sd,
                                  seq_sys, profile1, profile2,
                                  ins_prob1, del_prob1, ins_prob2, del_prob2,
                                  frag_len_min, frag_len_max,
                                  variant_probs, barcodes, pcr_dups) {

    # Check for improper argument types:
    check_ill_args(seq_object, n_reads, read_length, paired,
                   frag_mean, frag_sd, seq_sys, profile1, profile2,
                   ins_prob1, del_prob1, ins_prob2, del_prob2,
                   frag_len_min, frag_len_max, variant_probs, barcodes, pcr_dups)

    # Change mean and SD to shape and scale of Gamma distribution:
    frag_len_shape <- (frag_mean / frag_sd)^2
    frag_len_scale <- frag_sd^2 / frag_mean

    if (is.null(frag_len_min)) frag_len_min <- read_length
    # (Because it's coerced to an unsigned 32-bit integer, values >= 2^32 will
    #  make the integer "overflow" and change to something you don't want.)
    if (is.null(frag_len_max) || frag_len_max > 2^32 - 1) frag_len_max <- 2^32 - 1
    if (frag_len_min > frag_len_max) {
        stop("\nFragment length min can't be less than the max. ",
             "For computational reasons, it should also be < 2^32.", call. = FALSE)
    }
    if (is.null(variant_probs) && inherits(seq_object, "variants")) {
        variant_probs <- rep(1, seq_object$n_vars())
    }
    if (is.null(barcodes) && inherits(seq_object, "variants")) {
        barcodes <- rep("", seq_object$n_vars())
    } else if (is.null(barcodes)) barcodes <- ""

    sampler <- NULL
    if (paired) {
        prof_info1 <- read_profile(profile1, seq_sys, read_length, 1)
        prof_info2 <- read_profile(profile2, seq_sys, read_length, 2)
        if (inherits(seq_object, "ref_genome")) {
            sampler <- create_ref_ill_pe(seq_object$genome,
                                         frag_len_shape, frag_len_scale,
                                         frag_len_min, frag_len_max,
                                         prof_info1$qual_probs, prof_info1$quals,
                                         ins_prob1, del_prob1,
                                         prof_info2$qual_probs, prof_info2$quals,
                                         ins_prob2, del_prob2,
                                         barcodes)
        } else if (inherits(seq_object, "variants")) {
            sampler <- create_var_ill_pe(seq_object$genomes, variant_probs,
                                         frag_len_shape, frag_len_scale,
                                         frag_len_min, frag_len_max,
                                         prof_info1$qual_probs, prof_info1$quals,
                                         ins_prob1, del_prob1,
                                         prof_info2$qual_probs, prof_info2$quals,
                                         ins_prob2, del_prob2,
                                         barcodes)
        }
    } else {
        prof_info <- list(read_profile(profile1, seq_sys, read_length, 1))
        if (inherits(seq_object, "ref_genome")) {
            sampler <- create_ref_ill_se(seq_object$genome,
                                         frag_len_shape, frag_len_scale,
                                         frag_len_min, frag_len_max,
                                         prof_info$qual_probs, prof_info$quals,
                                         ins_prob1, del_prob1,
                                         barcodes)
        } else if (inherits(seq_object, "variants")) {
            sampler <- create_var_ill_se(seq_object$genomes, variant_probs,
                                         frag_len_shape, frag_len_scale,
                                         frag_len_min, frag_len_max,
                                         prof_info$qual_probs, prof_info$quals,
                                         ins_prob1, del_prob1,
                                         barcodes)
        }
    }

    return(sampler)
}




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
#' ARGUMENTS FOR FULL ILLUMINA SEQUENCING FUNCTION:
#'
#' seq_object,
#' out_prefixes,
#' n_reads,
#' read_length,
#' paired,
#' frag_mean,
#' frag_sd,
#' seq_sys = NULL,
#' profile1 = NULL,
#' profile2 = NULL,
#' ins_prob1 = 0.00009,
#' del_prob1 = 0.00011,
#' ins_prob2 = 0.00015,
#' del_prob2 = 0.00023,
#' frag_len_min = NULL,
#' frag_len_max = NULL,
#' variant_probs = NULL,
#' barcodes = NULL,
#' pcr_dups = 0.02,
#' id_info = list(),
#' n_cores
#'


#' Organize info for ID line in FASTQ file.
#'
#' @param instrument Instrument ID (string with the following characters allowed:
#'     a–z, A–Z, 0–9 and underscore). Defaults to `"SIM"`.
#' @param run_number Run number on instrument (numeric). Defaults to `1`.
#' @param flowcell_ID ID for flowcell (string with the following characters allowed:
#'     a–z, A–Z, 0–9). Defaults to `"FCX"`.
#' @param lane Lane number (numeric). Defaults to `1`.
#' @param tile Tile number (numeric). Defaults to `15`.
#' @param x_pos Run number on instrument (numeric). Defaults to `6329`.
#' @param y_pos X coordinate of cluster (numeric). Defaults to `1045`.
#' @param read Read number (numeric). Note that changing this doesn't do anything
#'     because it's set by the sequencing sampler as reads are created. Defaults to `1`.
#' @param is_filtered `"Y"` if the read is filtered (did not pass), `"N"` otherwise.
#'     Defaults to `"N"`.
#' @param control_number `0` when none of the control bits are on, otherwise it is
#'     an even number. On HiSeq X systems, control specification is not performed
#'     and this number is always `0`. Defaults to `0`.
#' @param sample_number Sample number from sample sheet (numeric). Defaults to `1`.
#'
#' @noRd
#'
#'
id_line_info <- function(...) {

    pars <- list(...)

    err_msg <- function(par, ...) {
        stop(sprintf("\nFor Illumina ID-line info, argument `%s` must be %s.",
                     par, paste(...)), call. = FALSE)
    }

    out <- list(
        instrument = "SIM",
        run_number = 1,
        flowcell_ID = "FCX",
        lane = 1,
        tile = 15,
        x_pos = 6329,
        y_pos = 1045,
        read = 1,
        is_filtered = "N",
        control_number = 0,
        sample_number = 1
    )

    # Check for nonsense if any parameters are provided, then insert them into output
    if (length(pars) > 0) {
        for (x in c("run_number", "lane", "tile", "x_pos", "y_pos", "read",
                    "sample_number")) {
            if (!is.null(pars[[x]]) && !single_integer(pars[[x]], 0)) {
                err_msg(x, "a single integer >= 0")
            }
        }
        for (x in c("instrument", "flowcell_ID")) {
            if (!is.null(pars[[x]])) {
                if (!is_type(pars[[x]], "character", 1)) err_msg(x, "a single string")
                allowed <- unlist(strsplit(c(letters, LETTERS, paste(0:9)), ""))
                if (x == "instrument") {
                    allowed <- c(allowed, "_")
                    allowed_msg <- "a–z, A–Z, 0–9 and underscore"
                } else allowed_msg <- "a–z, A–Z, and 0–9"
                if (!all(strsplit(pars[[x]], "")[[1]] %in% allowed)) {
                    err_msg(x, "a single string containing only", allowed_msg)
                }
            }
        }
        if (!is.null(pars$is_filtered)) {
            if (!is_type(pars$is_filtered, "character", 1) ||
                !pars$is_filtered %in% c("Y", "N")) {
                err_msg("is_filtered", "a single string of \"Y\" or \"N\"")
            }
        }
        if (!is.null(pars$control_number)) {
            if (!single_integer(pars$control_number, 0) ||
                pars$control_number %% 2 != 0) {
                err_msg("control_number", "a single even integer >= 0.")
            }
        }
        # Now put into output:
        for (x in names(out)) {
            if (!is.null(pars[[x]])) out[[x]] <- pars[[x]]
        }
    }

    return(out)

}



