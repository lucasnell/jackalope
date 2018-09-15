

# Cleaning up when the package is unloaded
.onUnload <- function (libpath) {
    library.dynam.unload("gemino", libpath)
}

# Check for a single, whole number, perhaps in range
single_whole_number <- function(x, .min, .max) {
    if (is.null(x)) return(FALSE)
    bool <- is.numeric(x) & length(x) == 1 & x %% 1 == 0
    if (!missing(.min)) bool <- bool & x >= .min
    if (!missing(.max)) bool <- bool & x <= .max
    return(bool)
}
# Check for a single number, perhaps in range
single_number <- function(x, .min, .max) {
    if (is.null(x)) return(FALSE)
    bool <- is.numeric(x) & length(x) == 1
    if (!missing(.min)) bool <- bool & x >= .min
    if (!missing(.max)) bool <- bool & x <= .max
    return(bool)
}
