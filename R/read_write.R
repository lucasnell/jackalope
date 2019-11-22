
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
#' @param cut_names Boolean for whether to cut chromosome names at the first space.
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


#' Write a `ref_genome` or `variants` object to a FASTA file.
#'
#' This file produces 1 FASTA file for a `ref_genome` object and one file
#' for each variant in a `variants` object.
#'
#' @param obj A `ref_genome` or `variants` object.
#' @param out_prefix Prefix for the output file.
#' @param compress Logical specifying whether or not to compress output file, or
#'     an integer specifying the level of compression, from 1 to 9.
#'     If `TRUE`, a compression level of `6` is used.
#'     Defaults to `FALSE`.
#' @param comp_method Character specifying which type of compression to use if any
#'     is desired. Options include `"gzip"` and `"bgzip"`.
#'     This is ignored if `compress` is `FALSE`. Defaults to `"bgzip"`.
#' @param text_width The number of characters per line in the output fasta file.
#'     Defaults to `80`.
#' @param show_progress Logical for whether to show a progress bar.
#'     Defaults to `FALSE`.
#' @param n_threads Number of threads to use if writing from a `variants` object.
#'     Threads are split among variants, so it's not useful to provide more threads
#'     than variants.
#'     This argument is ignored if `obj` is a `ref_genome` object, or if
#'     OpenMP is not enabled.
#'     Defaults to `1`.
#' @param overwrite Logical for whether to overwrite existing file(s) of the
#'     same name, if they exist. Defaults to `FALSE`.
#'
#' @return `NULL`
#'
#' @export
#'
write_fasta <- function(obj, out_prefix,
                        compress = FALSE,
                        comp_method = "bgzip",
                        text_width = 80,
                        show_progress = FALSE,
                        n_threads = 1,
                        overwrite = FALSE) {

    if (!inherits(obj, c("ref_genome", "variants"))) {
        err_msg("write_fasta", "obj", "a \"ref_genome\" or \"variants\" object")
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
    if (inherits(obj, "ref_genome")) {
        if (!inherits(obj$ptr(), "externalptr")) {
            stop("\nThe `ptr` method in the `obj` argument supplied to ",
                 "`write_fasta` should return an external pointer when the `obj` ",
                 "argument is of class \"ref_genome\".",
                 call. = TRUE)
        }
        check_file_existence(paste0(out_prefix, ".fa"), compress, overwrite)
        invisible(write_ref_fasta(out_prefix, obj$ptr(), text_width,
                                  compress, comp_method, show_progress))
    } else {
        if (!inherits(obj$ptr(), "externalptr")) {
            stop("\nThe `ptr` method in the `obj` argument supplied to ",
                 "`write_fasta` should return an external pointer when the `obj` ",
                 "argument is of class \"variants\".",
                 call. = TRUE)
        }
        check_file_existence(paste0(out_prefix, "__", obj$var_names(), ".fa"),
                             compress, overwrite)
        invisible(write_vars_fasta(out_prefix, obj$ptr(), text_width,
                                   compress, comp_method, n_threads, show_progress))
    }
    return(invisible(NULL))
}








# VCF ----



#' Write variant info from a \code{variants} object to a VCF file.
#'
#' Compression in this function always uses `"bgzip"` for compatibility with `"tabix"`.
#'
#' @param vars A \code{variants} object.
#' @inheritParams write_fasta
#' @param sample_matrix Matrix to specify how haploid variants are grouped into samples
#'     if samples are not haploid. There should be one row for each sample, and
#'     each row should contain indices or names for the variants present in that sample.
#'     Indices/names for variants cannot be repeated.
#'     Instead of repeating indices here, you should use the `dup_vars`
#'     method of the `variants` class to duplicate the necessary variant(s).
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
    if (!inherits(vars$ptr(), "externalptr")) {
        stop("\nThe `ptr` method in the `vars` argument supplied to ",
             "`write_vcf` should return an external pointer.",
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
    # Check for duplicates in `sample_matrix`:
    if (any(duplicated(as.numeric(sample_matrix)))) {
        stop("\nThe `sample_matrix` argument to the `write_vcf` function contained ",
             "duplicates.", call. = FALSE)
    }
    if (!is_type(show_progress, "logical", 1)) {
        err_msg("write_fasta", "show_progress", "a single logical")
    }
    if (!is_type(overwrite, "logical", 1)) {
        err_msg("write_fasta", "overwrite", "a single logical")
    }

    check_file_existence(paste0(out_prefix, ".vcf"), compress, overwrite)

    write_vcf_cpp(out_prefix, compress, vars$ptr(), sample_matrix, show_progress)

    return(invisible(NULL))
}

