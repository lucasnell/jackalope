

# Cleaning up when the package is unloaded
.onUnload <- function (libpath) {
    library.dynam.unload("gemino", libpath)
}
