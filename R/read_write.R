

#' Read a fasta file to an external pointer object.
#'
#' Accepts uncompressed and gzipped fasta files.
#'
#' @param fasta_file File name of the fasta file.
#' @param fai_file File name of the fasta file.
#'     Providing this argument speeds up the reading process significantly.
#'     Defaults to \code{NULL}, which indicates the fasta file is not indexed.
#' @param cut_names Boolean for whether to cut sequence names at the first space.
#'     Defaults to \code{TRUE}.
#' @param rm_soft_mask Boolean for whether to remove soft-masking by making
#'    sequences all uppercase. Defaults to \code{TRUE}.
#'
#' @return An external pointer to a \code{RefGenome} object in C++.
#'
#' @export
#'
#'
read_fasta <- function(fasta_file, fai_file = NULL,
                       cut_names = TRUE, rm_soft_mask = TRUE) {

    if (is.null(fai_file)) {
        ptr <- read_fasta_noind(fasta_file, cut_names, rm_soft_mask)
    } else {
        ptr <- read_fasta_ind(fasta_file, fai_file, rm_soft_mask)
    }

    return(ptr)
}


#' Write to FASTA file.
#'
#' @param file_name File name of the output fasta file.
#' @param ptr External pointer to a \code{RefGenome} C++ object.
#' @param text_width The number of characters per line in the output fasta file.
#'     Defaults to 80.
#' @param compression Type of compression. Takes either \code{"none"} or \code{"gzip"}.
#'     Defaults to \code{"none"}.
#'
#' @return Nothing.
#'
#' @export
#'
write_fasta <- function(file_name, ptr, text_width = 80, compression = "none") {
    if (compression == "gzip" || compression == "gz") {
        invisible(write_fasta_gz(file_name, ptr, text_width))
    } else if (compression == "none") {
        invisible(write_fasta_fa(file_name, ptr, text_width))
    } else {
        stop("Unknown compression type. Only 'none' or 'gzip' are implemented.");
    }
    return(invisible(NULL))
}
