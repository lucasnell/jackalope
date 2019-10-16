

# Cleaning up when the package is unloaded
.onUnload <- function (libpath) {
    library.dynam.unload("jackalope", libpath)
}

# For avoiding warnings for comparing nonsensible inputs
comparable <- function(x) {
    return(!is.null(x) && (is.vector(x) || is.matrix(x)) && !any(is.na(x)))
}

# Check for a single, whole number, perhaps in range
single_integer <- function(x, .min, .max) {
    if (!comparable(x)) return(FALSE)
    bool <- is.numeric(x) && length(x) == 1 && x %% 1 == 0
    if (!missing(.min)) bool <- bool && x >= .min
    if (!missing(.max)) bool <- bool && x <= .max
    return(bool)
}
# Check for a single number, perhaps in range
single_number <- function(x, .min, .max) {
    if (!comparable(x)) return(FALSE)
    bool <- is.numeric(x) && length(x) == 1
    if (!missing(.min)) bool <- bool && x >= .min
    if (!missing(.max)) bool <- bool && x <= .max
    return(bool)
}
# Check for a vector of positive numbers that sums to > 0 (used often for relative rates)
positive_vector <- function(x, zero_comp = `>=`) {
    if (!comparable(x)) return(FALSE)
    bool <- inherits(x, c("integer", "numeric")) && zero_comp(sum(x), 0) && all(x >= 0)
    return(bool)
}
is_type <- function(x, type, L = NULL) {
    if (!comparable(x)) return(FALSE)
    if (!inherits(x, type)) return(FALSE)
    if (!is.null(L) && !length(x) %in% L) return(FALSE)
    return(TRUE)
}

#' Standard way to show error messages (also to make input-checking less verbose).
#'
#' @noRd
#'
err_msg <- function(fxn, par, ...) {
    stop(sprintf("\nFor the `%s` function in jackalope, argument `%s` must be %s.",
                 fxn, par, paste(...)), call. = FALSE)
}




#' Check for whether file(s) already exist, return error depending on `overwrite` arg.
#'
#'
#' @noRd
#'
check_file_existence <- function(file_names, compress, overwrite) {

    file_names <- path.expand(file_names)

    if (compress) file_names <- paste0(file_names, ".gz")

    dir_names <- unique(dirname(file_names))

    for (d in dir_names) {
        if (!dir.exists(d)) dir.create(d, recursive = TRUE)
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
