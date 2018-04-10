#
# Digest genome(s)
#





# ----------------
# Digest a variants object
# ----------------

digest_variants <- function(variants_obj, enzyme_names, n_cores,
                            chunk_size, enz_list, in_place) {


    enzyme_sites <- c(enz_list[enzyme_names], recursive = TRUE, use.names = FALSE)
    len5s <- sapply(seq(1, length(enzyme_sites), 2), function(i) nchar(enzyme_sites[i]))
    bind_sites <- sapply(seq(1, length(enzyme_sites), 2),
                         function(i) cpp_merge_str(enzyme_sites[i:(i+1)]))

    indiv_lists <- digest_all_variants_all_scaffs(
        variants_obj$reference, variants_obj$variant_set,
        bind_sites, len5s, chunk_size, n_cores)

    if (in_place) {
        variants_obj$digests <- indiv_lists
        variants_obj$binding_sites <- enzyme_sites
        invisible(NULL)
    } else {
        variants_out <- variants_obj$clone()
        variants_out$digests <- indiv_lists
        variants_out$binding_sites <- enzyme_sites
        return(variants_out)
    }
}











# ----------------
# Digest a dna_set object
# ----------------

digest_reference <- function(dna_set_in, enzyme_names, n_cores, enz_list, in_place) {

    enzyme_sites <- c(enz_list[enzyme_names], recursive = TRUE, use.names = FALSE)
    len5s <- sapply(seq(1, length(enzyme_sites), 2), function(i) nchar(enzyme_sites[i]))
    bind_sites <- sapply(seq(1, length(enzyme_sites), 2),
                         function(i) cpp_merge_str(enzyme_sites[i:(i+1)]))

    dig_list <- digest_ref_scaffs(dna_set_in$sequence_set, bind_sites, len5s, n_cores)

    if (in_place) {
        dna_set_in$digests <- dig_list
        dna_set_in$binding_sites <- enzyme_sites
        invisible(NULL)
    } else {
        dna_set_out <- dna_set_in$clone()
        dna_set_out$digests <- dig_list
        dna_set_out$binding_sites <- enzyme_sites
        return(dna_set_out)
    }
}




#'
#' Digest genome(s).
#'
#' \emph{Note:} This will override any digestions currently in place in the
#' object. If you want to add a new digestion, re-run this function with the names
#' of all enzymes you're interested in included in the \code{enzyme_names} argument.
#'
#'
#' @param object Either a \code{dna_set} or \code{variants} object.
#' @param enzyme_names Name of enzyme(s).
#' @param n_cores Number of cores to use for parallel processing. This argument is
#'     ignored if OpenMP is not enabled. Defaults to \code{1}.
#' @param chunk_size The size of chunks to break up scaffolds into when digesting
#'     a \code{variants} object.
#'     (This argument is ignored if digesting a \code{dna_set}.)
#'     Changing this might affect performance, for better or worse.
#'     The default worked best on my computer. Defaults to \code{1000}.
#' @param enz_list List of enzymes with binding sites. Default is the internal
#'     \code{binding_sites} list (see \code{\link{binding_sites}}).
#' @param in_place Boolean for whether to edit the object in place without
#'     making a new copy. Defaults to \code{FALSE}.
#'
#' @return If \code{in_place == FALSE}, a \code{variants} or \code{dna_set} object
#'     with the \code{digests} field filled in.
#'     If \code{in_place == TRUE}, it returns \code{NULL}, but it changes the input
#'     object in place.
#'
#'
#' @export
#'
#' @examples
#'
#' ref_genome <- dna_set$new(rando_seqs(100, mean_len = 1e3, sd_len = 1e2))
#' digest(ref_genome, 'ApeKI', n_cores = 1, in_place = TRUE)
#' ref_genome
#'
#' variants_obj <- make_variants(ref_genome, n_vars = 10)
#' variants_obj
#'
#' # Returns a new variants object
#' digest(variants_obj, 'ApeKI')
#'
#' # Returns nothing, but changes variants_obj object
#' digest(variants_obj, 'AscI', in_place = TRUE)
#' # To see the changes...
#' variants_obj
#'
digest <- function(object, enzyme_names, n_cores = 1,
                   chunk_size = 1000,
                   enz_list = binding_sites, in_place = FALSE) {

    if (missing(enzyme_names)) stop("enzyme_names required")
    if (missing(object)) stop("object required")

    if (any(!enzyme_names %in% names(enz_list))) {
        cat(paste(names(enz_list), collapse = "  "))
        stop(paste("One or more of the enzyme_names input is not present in the",
                   "names of the enz_list argument.",
                   "See output above for the names that are present."))
    }

    if (is(object, 'dna_set')) {
        digest_reference(object, enzyme_names, n_cores, enz_list, in_place)
    } else if (is(object, 'variants')) {
        digest_variants(object, enzyme_names, n_cores, chunk_size, enz_list, in_place)
    } else {
        stop("Input object must be a dna_set or variants object.")
    }
    invisible(NULL)
}
