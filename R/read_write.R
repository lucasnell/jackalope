

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

    ref_obj <- ref_genome$new(ptr)

    return(ref_obj)
}


#' Write a \code{ref_genome} object to a FASTA file.
#'
#' @param ref_obj A \code{ref_genome} object.
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
write_fasta <- function(ref_obj, file_name, text_width = 80, compress = FALSE) {

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
    if (!inherits(ref_obj, "ref_genome")) {
        stop("\nThe ref_obj argument supplied to write_fasta is not a ref_genome object.",
             call. = FALSE)
    }
    if (!inherits(ref_obj$genome, "externalptr")) {
        stop("\nThe genome field in the ref_genome object supplied to write_fasta ",
             "is not an external pointer.",
             call. = TRUE)
    }
    if (compress) {
        invisible(write_fasta_gz(file_name, ref_obj$genome, text_width))
    } else {
        invisible(write_fasta_fa(file_name, ref_obj$genome, text_width))
    }
    return(invisible(NULL))
}
