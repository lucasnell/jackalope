

# Cleaning up when the package is unloaded
.onUnload <- function (libpath) {
    library.dynam.unload("jackal", libpath)
}

# Check for a single, whole number, perhaps in range
single_integer <- function(x, .min, .max) {
    if (is.null(x) || any(is.na(x))) return(FALSE)
    bool <- is.numeric(x) & length(x) == 1 & x %% 1 == 0
    if (!missing(.min)) bool <- bool & x >= .min
    if (!missing(.max)) bool <- bool & x <= .max
    return(bool)
}
# Check for a single number, perhaps in range
single_number <- function(x, .min, .max) {
    if (is.null(x) || any(is.na(x))) return(FALSE)
    bool <- is.numeric(x) & length(x) == 1
    if (!missing(.min)) bool <- bool & x >= .min
    if (!missing(.max)) bool <- bool & x <= .max
    return(bool)
}
is_type <- function(x, type, L = NULL) {
    if (is.null(x) || any(is.na(x))) return(FALSE)
    if (!inherits(x, type)) return(FALSE)
    if (!is.null(L) && !length(x) %in% L) return(FALSE)
    return(TRUE)
}

#' Standard way to show error messages (also to make input-checking less verbose).
#'
#' @noRd
#'
err_msg <- function(fxn, par, ...) {
    stop(sprintf("\nFor the `%s` function in jackal, argument `%s` must be %s.",
                 fxn, par, paste(...)), call. = FALSE)
}
