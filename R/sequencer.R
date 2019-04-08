

#' Check for whether a fastq file already exists.
#'
#'
#' @noRd
#'
check_fastq <- function(out_prefix, compress, overwrite, n_read_ends = 1) {

    out_prefix <- path.expand(out_prefix)

    file_names <- paste0(out_prefix, "_R", 1:n_read_ends, ".fq")
    if (compress) file_names <- paste0(file_names, ".gz")

    dir_names <- unique(dirname(file_names))

    for (d in dir_names) {
        if (!dir.exists(d)) dir.create(d)
    }

    if (!overwrite) {
        for (f in file_names) {
            if (file.exists(f)) {
                stop("\nFile ", paste(f), " already exists.", call. = FALSE)
            }
        }

    }

    invisible(NULL)

}



