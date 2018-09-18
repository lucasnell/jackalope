
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
    if (single_string(method_info)) {

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
             "(using `$change_names()` method) to have the same names as the VCF file. ",
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




