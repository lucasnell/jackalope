
# FASTA ----


#' Read a fasta file.
#'
#' Accepts uncompressed and gzipped fasta files.
#'
#' @param fasta_file File name of the fasta file.
#' @param fai_file File name of the fasta index file.
#'     Providing this argument speeds up the reading process significantly.
#'     Defaults to \code{NULL}, which indicates the fasta file is not indexed.
#' @param cut_names Boolean for whether to cut sequence names at the first space.
#'     This argument is ignored if \code{fai_file} is not \code{NULL}.
#'     Defaults to \code{TRUE}.
#' @param rm_soft_mask Boolean for whether to remove soft-masking by making
#'    sequences all uppercase. Defaults to \code{TRUE}.
#'
#' @return A \code{\link{ref_genome}} object.
#'
#' @export
#'
#'
read_fasta <- function(fasta_file, fai_file = NULL,
                       cut_names = TRUE, rm_soft_mask = TRUE) {


    if (!inherits(cut_names, "logical") | length(cut_names) != 1) {
        stop("\nThe cut_names argument supplied to read_fasta is not a logical",
             "vector of length one.",
             call. = FALSE)
    }
    if (!inherits(rm_soft_mask, "logical") | length(rm_soft_mask) != 1) {
        stop("\nThe rm_soft_mask argument supplied to read_fasta is not a logical",
             "vector of length one.",
             call. = FALSE)
    }
    if (!inherits(fasta_file, "character") | length(fasta_file) != 1) {
        stop("\nThe fasta_file argument supplied to read_fasta is not a character ",
             "vector of length one.",
             call. = FALSE)
    }

    if (is.null(fai_file)) {
        ptr <- read_fasta_noind(fasta_file, cut_names, rm_soft_mask)
    } else {
        if (!inherits(fai_file, "character") | length(fai_file) != 1) {
            stop("\nThe fai_file argument supplied to read_fasta is neither NULL nor ",
                 "a character vector of length one ",
                 call. = FALSE)
        }
        ptr <- read_fasta_ind(fasta_file, fai_file, rm_soft_mask)
    }

    reference <- ref_genome$new(ptr)

    return(reference)
}


#' Write a \code{ref_genome} object to a FASTA file.
#'
#' @param reference A \code{ref_genome} object.
#' @param file_name File name of the output fasta file.
#' @param text_width The number of characters per line in the output fasta file.
#'     Defaults to \code{80}.
#' @param compress Boolean for whether to compress using \code{"gzip"}.
#'     Defaults to \code{FALSE}.
#'
#' @return \code{NULL}
#'
#' @export
#'
write_fasta <- function(reference, file_name, text_width = 80, compress = FALSE) {

    if (!inherits(compress, "logical") | length(compress) != 1) {
        stop("\nThe compress argument supplied to write_fasta is not a logical",
             "vector of length one.",
             call. = FALSE)
    }
    if (!is.numeric(text_width) | length(text_width) != 1) {
        stop("\nThe text_width argument supplied to write_fasta is not a numeric ",
             "vector of length one.",
             call. = FALSE)
    }
    if (!inherits(file_name, "character") | length(file_name) != 1) {
        stop("\nThe file_name argument supplied to write_fasta is not a character ",
             "vector of length one.",
             call. = FALSE)
    }
    if (!inherits(reference, "ref_genome")) {
        stop("\nThe reference argument supplied to write_fasta is not a ref_genome ",
             "object.",
             call. = FALSE)
    }
    if (!inherits(reference$genome, "externalptr")) {
        stop("\nThe genome field in the ref_genome object supplied to write_fasta ",
             "is not an external pointer.",
             call. = TRUE)
    }
    if (compress) {
        invisible(write_fasta_gz(file_name, reference$genome, text_width))
    } else {
        invisible(write_fasta_fa(file_name, reference$genome, text_width))
    }
    return(invisible(NULL))
}








# VCF ----


#' Read VCF file using package `vcfR`.
#'
#' This function is only to be used inside `create_variants`.
#'
#' @inheritParams write_fasta
#' @param method_info List or string passed from `create_variants`. If a string, it
#'     should be the file name of VCF file.
#'     If a list, it must contain the VCF file name in the `file` field.
#'     The list can also optionally contain
#'     a boolean for whether to print chromosomes present in the VCF file
#'     (`print_names`),
#'     and other parameters passed to \code{\link[vcfR]{read.vcfR}}.
#'
#' @noRd
#'
read_vcf <- function(reference, method_info) {

    if (!requireNamespace("vcfR", quietly = TRUE)) {
        stop("\nPackage \"vcfR\" is needed for reading VCF files. ",
             "Please install it.",
             call. = FALSE)
    }

    err <- FALSE
    if (is_type(method_info, "character", 1)) {

        print_names <- FALSE
        read_args <- list(file = method_info, verbose = FALSE)

    } else if (inherits(method_info, "list")) {

        if (is.null(method_info$file)) {
            err <- TRUE
        } else {

            if (is.null(method_info$print_names)) {
                print_names <- FALSE
            } else print_names <- method_info$print_names

            read_args <- method_info
            # No longer needed in list (bc not used in `vcfR::readvcfR`):
            read_args$print_names <- NULL
            if (is.null(read_args$verbose)) read_args$verbose <- FALSE
            if (!all(names(read_args) %in% names(formals(vcfR::read.vcfR)))) {
                err <- TRUE
            }
        }

    } else err <- TRUE

    if (err) {
        stop("\nIf method = \"vcf\" in `create_variants`, ",
             "the `method_info` arg must be a single string specifying the ",
             "filename for the VCF file, or a list with the vcf file and arguments ",
             "for `vcfR::read.vcfR`. ",
             "See \"Method arguments\" section under `?create_variants`.",
             call. = FALSE)
    }

    vcf <- do.call(vcfR::read.vcfR, read_args)

    chrom <- vcf@fix[,"CHROM"]
    pos <- as.integer(vcf@fix[,"POS"]) - 1  # -1 is to convert to C++ indices
    ref_seq <- vcf@fix[,"REF"]

    if (length(chrom) != length(pos) | length(chrom) != length(ref_seq)) {
        stop("\nVCF not parsing correctly. ",
             "Vectors of chromosomes, positions, and reference-sequences aren't ",
             "all the same length.",
             call. = FALSE)
    }

    seq_names <- view_ref_genome_seq_names(reference$genome)
    unq_chrom <- unique(chrom)

    if (!all(unq_chrom %in% seq_names)) {
        if (print_names) print(unq_chrom)
        stop("\nSequence name(s) in VCF file don't match those in the ",
             "`ref_genome` object. ",
             "It's probably easiest to manually change the `ref_genome` object ",
             "(using `$set_names()` method) to have the same names as the VCF file. ",
             "Re-run this function with `print_names = TRUE` to see the VCF-file names.",
             call. = FALSE)
    }

    # Converts items in `chrom` to 0-based indices of sequences in ref. genome
    chrom_inds <- match(chrom, seq_names) - 1

    haps <- vcfR::extract.haps(vcf, unphased_as_NA = FALSE, verbose = FALSE)

    var_names <- colnames(haps)
    colnames(haps) <- NULL
    rownames(haps) <- NULL

    # So I don't have to deal with NAs:
    haps[is.na(haps)] <- ""

    if (nrow(haps) != length(pos)) {
        stop("\nVCF not parsing correctly. ",
             "Number of haplotypes doesn't match with number of positions.",
             call. = FALSE)
    }

    # Split into list for easier processing in `read_vcfr`
    haps <- split(haps, row(haps))

    variants_ptr <- read_vcfr(reference$genome, var_names,
                              haps, chrom_inds, pos, ref_seq)

    return(variants_ptr)

}






# Illumina profile ----


#' Find a package-supplied ART-style Illumina profile file.
#'
#' @inheritParams read_profile
#'
#' @return The filename to be read from `read_profile`.
#'
#' @noRd
#'
find_profile_file <- function(profile_name, read_length, read) {

    abbrevs <- list(GA1 = "Genome Analyzer I",
                    GA2 = "Genome Analyzer II",
                    HS10 = "HiSeq 1000",
                    HS20 = "HiSeq 2000",
                    HS25 = "HiSeq 2500",
                    MSv1 = "MiSeq v1",
                    MSv3 = "MiSeq v3",
                    MinS = "MiniSeq TruSeq",
                    NS50 = "NextSeq 500 v2",
                    HSXn = "HiSeqX v2.5 PCR free",
                    HSXt = "HiSeqX v2.5 TruSeq")
    if (profile_name %in% names(abbrevs)) {
        profile_name <- abbrevs[[profile_name]]
    }
    profile_lookup <- read.csv(text =
                                   "name,read_length,read,file_name
                               Genome Analyzer I,36,1,QUAL_DIST_ONE36
                               Genome Analyzer I,36,2,QUAL_DIST_TWO36
                               Genome Analyzer I,44,1,QUAL_DIST_ONE44
                               Genome Analyzer I,44,2,QUAL_DIST_TWO44
                               Genome Analyzer II,50,1,QUAL_DIST_ONE50
                               Genome Analyzer II,50,2,QUAL_DIST_TWO50
                               Genome Analyzer II,75,1,QUAL_DIST_ONE75
                               Genome Analyzer II,75,2,QUAL_DIST_TWO75
                               HiSeq 1000,100,1,QUAL_DIST_ONE100
                               HiSeq 1000,100,2,QUAL_DIST_TWO100
                               MiSeq v1,250,1,QUAL_DIST_ONE250
                               MiSeq v1,250,2,QUAL_DIST_TWO250
                               HiSeq 2000,100,1,HiSeq2000L100R1
                               HiSeq 2000,100,2,HiSeq2000L100R2
                               HiSeq 2500,125,1,HiSeq2500L125R1
                               HiSeq 2500,125,2,HiSeq2500L125R2
                               HiSeq 2500,150,1,HiSeq2500L150R1
                               HiSeq 2500,150,2,HiSeq2500L150R2
                               MiniSeq TruSeq,50,1,MiniSeqTruSeq_L50
                               MiSeq v3,250,1,MiSeqv3_L250_R1
                               MiSeq v3,250,2,MiSeqv3_L250_R2
                               NextSeq 500 v2,75,1,NextSeq500v2_L75_R1
                               NextSeq 500 v2,75,2,NextSeq500v2_L75_R2
                               HiSeqX v2.5 PCR free,150,1,HiSeqXPCRfree_L150_R1
                               HiSeqX v2.5 PCR free,150,2,HiSeqXPCRfree_L150_R2
                               HiSeqX v2.5 TruSeq,150,1,HiSeqXtruSeq_L150_R1
                               HiSeqX v2.5 TruSeq,150,2,HiSeqXtruSeq_L150_R2",
                     stringsAsFactors = FALSE)
        profile_lookup$name <- trimws(profile_lookup$name)
        profile_lookup <- profile_lookup[order(profile_lookup$read_length),]
        criteria <- profile_lookup$name == profile_name
        if (sum(criteria) == 0) {
            stop("\nThe desired Illumina platform name isn't available.", call. = FALSE)
        }
        criteria <- criteria & profile_lookup$read_length >= read_length
        if (sum(criteria) == 0) {
            stop("\nFor the desired Illumina platform, this package doesn't have ",
                 "a read length that's as long as you want.", call. = FALSE)
        }
        criteria <- criteria & profile_lookup$read == read
        if (sum(criteria) == 0) {
            stop("\nFor the desired Illumina platform of read length ",
                 paste(read_length), ", this package only has ",
                 "a profile for read number ",
                 paste(ifelse(read == 1, 2, 1)), ".", call. = FALSE)
        }
        profile_name <- paste0(profile_lookup$file_name[criteria][1], ".txt")
        profile_name <- system.file("art_profiles", profile_name, package = "gemino",
                                    mustWork = TRUE)
}






#' Read an ART-style Illumina profile file.
#'
#' @param profile_name File name of ART-style profile file for Illumina reads.
#'     If this argument doesn't end in `.txt`, then it's assumed it's the name of
#'     an Illumina platform with profile(s) saved in the package.
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
read_profile <- function(profile_name, read_length, read) {

    if (!grepl(".txt$", profile_name)) {
        profile_name <- find_profile_file(profile_name, read_length, read)
    }

    profile_str <- readLines(profile_name)
    profile_str <- profile_str[profile_str != ""]
    profile_str <- profile_str[!grepl("^\\.|^#|^N", profile_str)]
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


    pos <- sapply(profile_info, function(x) x$pos)
    if (min(pos) > 0) stop("\nMinimum profile position should be zero.", call. = FALSE)
    if (max(pos) < (read_length-1)) {
        stop("\nMaximum profile position should be >= read_length - 1.", call. = FALSE)
    }

    nts <- sapply(profile_info, function(x) x$nt)

    mis_probs <- rep(list(NA), 4)
    quals <- rep(list(NA), 4)
    names(mis_probs) <- names(quals) <- c("T", "C", "A", "G")

    for (nt in names(mis_probs)) {
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
        # Remove unnecessary positions if necessary:
        if (length(probs_nt) > read_length) {
            probs_nt <- probs_nt[1:read_length]
            qual_nt <- qual_nt[1:read_length]
        }
        mis_probs[[nt]] <- probs_nt
        quals[[nt]] <- qual_nt
    }

    # Names no longer needed
    names(mis_probs) <- names(quals) <- NULL

    return(list(mis_probs = mis_probs, quals = quals))

}

