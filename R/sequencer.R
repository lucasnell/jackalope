
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



#' Check for whether a fastq file already exists.
#'
#'
#' @noRd
#'
check_fastq <- function(out_prefix, n_read_ends = 1) {

    out_prefix <- path.expand(out_prefix)

    file_names <- paste0(out_prefix, "_R", 1:n_read_ends, ".fq.gz")
    dir_names <- unique(dirname(file_names))

    for (d in dir_names) {
        if (!dir.exists(d)) dir.create(d)
    }

    for (f in file_names) {
        if (file.exists(f)) stop("\nFile ", paste(f), " already exists.", call. = FALSE)
    }

    invisible(NULL)

}



