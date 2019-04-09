
# FASTA ----


#' Read a fasta file.
#'
#' Accepts uncompressed and gzipped fasta files.
#'
#' @param fasta_files File name(s) of the fasta file(s).
#' @param fai_files File name(s) of the fasta index file(s).
#'     Providing this argument speeds up the reading process significantly.
#'     If this argument is provided, it must be the same length as the `fasta_files`
#'     argument.
#'     Defaults to \code{NULL}, which indicates the fasta file(s) is/are not indexed.
#' @param cut_names Boolean for whether to cut sequence names at the first space.
#'     This argument is ignored if \code{fai_file} is not \code{NULL}.
#'     Defaults to \code{FALSE}.
#'
#' @return A \code{\link{ref_genome}} object.
#'
#' @export
#'
#'
read_fasta <- function(fasta_files, fai_files = NULL,
                       cut_names = FALSE) {


    if (!is_type(fasta_files, "character")) {
        err_msg("read_fasta", "fasta_files", "a character vector")
    }
    if (!is.null(fai_files) && !is_type(fai_files, "character", length(fasta_files))) {
        err_msg("read_fasta", "fai_files", "NULL or a character vector of the same",
                "length as the `fasta_files` argument")
    }
    if (!is_type(cut_names, "logical", 1)) {
        err_msg("read_fasta", "cut_names", "a single logical")
    }

    # For now I'm forcing the users to remove soft-masking
    rm_soft_mask <- TRUE

    if (is.null(fai_files)) {
        ptr <- read_fasta_noind(fasta_files, cut_names, rm_soft_mask)
    } else {
        ptr <- read_fasta_ind(fasta_files, fai_files, rm_soft_mask)
    }

    reference <- ref_genome$new(ptr)

    return(reference)
}


#' Write a sequence object to a FASTA file.
#'
#' This file produces 1 FASTA file for a `ref_genome` object and one file
#' for each variant in a `variants` object.
#'
#' @param seq_obj A `ref_genome` or `variants` object.
#' @param out_prefix Prefix for the output file.
#' @param compress Logical specifying whether or not to compress output file, or
#'     an integer specifying the level of compression, from 1 to 9.
#'     If `TRUE`, a compression level of `6` is used.
#'     Defaults to `FALSE`.
#' @param comp_method Character specifying which type of compression to use if any
#'     is desired. Options include `"gzip"` and `"bgzip"`. Defaults to `"gzip"`.
#' @param text_width The number of characters per line in the output fasta file.
#'     Defaults to `80`.
#' @param show_progress Logical for whether to show a progress bar.
#'     Defaults to `FALSE`.
#' @param n_threads Number of threads to use if writing from a `variants` object.
#'     Threads are split among variants, so it's not useful to provide more threads
#'     than variants. This argument is ignored if `seq_obj` is a `ref_genome` object.
#'     Defaults to `1`.
#' @param overwrite Logical for whether to overwrite existing file(s) of the
#'     same name, if they exist. Defaults to `FALSE`.
#'
#' @return `NULL`
#'
#' @export
#'
write_fasta <- function(seq_obj, out_prefix,
                        compress = FALSE,
                        comp_method = "gzip",
                        text_width = 80,
                        show_progress = FALSE,
                        n_threads = 1,
                        overwrite = FALSE) {

    if (!inherits(seq_obj, c("ref_genome", "variants"))) {
        err_msg("write_fasta", "seq_obj", "a \"ref_genome\" or \"variants\" object")
    }
    if (!is_type(out_prefix, "character", 1)) {
        err_msg("write_fasta", "out_prefix", "a single string")
    }
    if (!is_type(compress, "logical", 1) && !single_integer(compress, 1, 9)) {
        err_msg("write_fasta", "compress", "a single logical or integer from 1 to 9")
    }
    if (is_type(compress, "logical", 1) && compress) compress <- 6 # default compression
    if (is_type(compress, "logical", 1) && !compress) compress <- 0 # no compression
    if (!is_type(comp_method, "character", 1) || !comp_method %in% c("gzip", "bgzip")) {
        err_msg("write_fasta", "comp_method", "\"gzip\" or \"bgzip\"")
    }
    if (!single_integer(text_width, 1)) {
        err_msg("write_fasta", "text_width", "a single integer >= 1")
    }
    if (!is_type(show_progress, "logical", 1)) {
        err_msg("write_fasta", "show_progress", "a single logical")
    }
    if (!single_integer(n_threads, 1)) {
        err_msg("write_fasta", "n_threads", "a single integer >= 1")
    }
    if (!is_type(overwrite, "logical", 1)) {
        err_msg("write_fasta", "overwrite", "a single logical")
    }
    if (inherits(seq_obj, "ref_genome")) {
        if (!inherits(seq_obj$genome, "externalptr")) {
            stop("\nThe `genome` field in the `seq_obj` argument supplied to ",
                 "`write_fasta` should be an external pointer when the `seq_obj` ",
                 "argument is of class \"ref_genome\".",
                 call. = TRUE)
        }
        invisible(write_ref_fasta(out_prefix, seq_obj$genome, text_width,
                                  compress, comp_method, show_progress))
    } else {
        if (!inherits(seq_obj$genomes, "externalptr")) {
            stop("\nThe `genomes` field in the `seq_obj` argument supplied to ",
                 "`write_fasta` should be an external pointer when the `seq_obj` ",
                 "argument is of class \"variants\".",
                 call. = TRUE)
        }
        invisible(write_vars_fasta(out_prefix, seq_obj$genomes, text_width,
                                   compress, comp_method, n_threads, show_progress))
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
    alts <- strsplit(vcf@fix[,"ALT"], ",")

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

    if (nrow(haps) != length(pos)) {
        stop("\nVCF not parsing correctly. ",
             "Number of haplotypes doesn't match with number of positions.",
             call. = FALSE)
    }

    # I'm assuming NAs mean no mutation
    haps[is.na(haps)] <- ""

    # Split into list for easier processing in `read_vcfr`
    haps <- split(haps, row(haps))

    # We treat things differently if vcfR has output numbers rather than nucleotides.
    # (The below line should be TRUE when it outputs numbers.)
    if (any(!is.na(suppressWarnings(as.integer(do.call(c, haps)))))) {
        # Change string integers to actual genotypes:
        haps <-
            lapply(1:length(haps),
                   function(i) {
                       as.character(sapply(haps[[i]],
                                           function(j) {
                                               ifelse(j == "" | j == "0", "",
                                                      alts[[i]][as.integer(j)])
                                           }))
                   })
    }


    variants_ptr <- read_vcfr(reference$genome, var_names,
                              haps, chrom_inds, pos, ref_seq)

    return(variants_ptr)

}



#' Write variant info from a \code{variants} object to a VCF file.
#'
#' Compression in this function always uses `"bgzip"` for compatibility with `"tabix"`.
#'
#' @param vars A \code{variants} object.
#' @inheritParams write_fasta
#' @param sample_matrix Matrix to specify how haploid variants are grouped into samples
#'     if samples are not haploid. There should be one row for each sample, and
#'     each row should contain indices or names for the variants present in that sample.
#'     Indices/names for variants can be repeated across and within rows.
#'     The number of columns indicates the ploidy level: 2 columns for diploid,
#'     3 for triploid, 4 for tetraploid, and so on;
#'     there is no limit to the ploidy level.
#'     If this argument is `NULL`, it's assumed that each variant is its own
#'     separate sample.
#'     Defaults to `NULL`.
#'
#' @return \code{NULL}
#'
#' @export
#'
write_vcf <- function(vars,
                      out_prefix,
                      compress = FALSE,
                      sample_matrix = NULL,
                      show_progress = FALSE,
                      overwrite = FALSE) {

    if (!inherits(vars, "variants")) {
        err_msg("write_vcf", "vars", "a \"variants\" object")
    }
    if (!is_type(out_prefix, "character", 1)) {
        err_msg("write_vcf", "out_prefix", "a single string")
    }
    if (!is_type(compress, "logical", 1) && !single_integer(compress, 1, 9)) {
        err_msg("write_fasta", "compress", "a single logical or integer from 1 to 9")
    }
    if (is_type(compress, "logical", 1) && compress) compress <- 6 # default compression
    if (is_type(compress, "logical", 1) && !compress) compress <- 0 # no compression
    if (!inherits(vars$genomes, "externalptr")) {
        stop("\nThe `genomes` field in the `vars` argument supplied to ",
             "`write_vcf` should be an external pointer.",
             call. = TRUE)
    }
    if (!is.null(sample_matrix) && !inherits(sample_matrix, "matrix")) {
        err_msg("write_vcf", "sample_matrix", "NULL or a matrix")
    }
    if (is.null(sample_matrix)) {
        sample_matrix <- cbind(1:vars$n_vars())
    }
    # If they provided names rather than indices:
    if (inherits(sample_matrix[1,1], "character")) {
        var_names <- vars$var_names()
        not_found <- !unique(as.character(sample_matrix)) %in% var_names
        if (sum(not_found) > 0) {
            stop("\nThe `sample_matrix` argument to the `write_vcf` function had ",
                 "the following variant name(s) that weren't found in the ",
                 "input variants object: ",
                 paste(unique(as.character(sample_matrix))[not_found], collapse = ", "),
                 call. = FALSE)
        }
        sample_matrix <- apply(sample_matrix, 1:2, function(nm) which(var_names == nm))
    } else if (!inherits(sample_matrix[1,1], c("numeric", "integer"))) {
        err_msg("write_vcf", "sample_matrix", "NULL or a character, numeric, or",
                "integer matrix.")
    }
    if (!is_type(show_progress, "logical", 1)) {
        err_msg("write_fasta", "show_progress", "a single logical")
    }
    if (!is_type(overwrite, "logical", 1)) {
        err_msg("write_fasta", "overwrite", "a single logical")
    }


    write_vcf_cpp(out_prefix, compress, vars$genomes, sample_matrix, show_progress)

    return(invisible(NULL))
}

