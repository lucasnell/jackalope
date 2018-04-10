
#'
#' Filter cut sites based on resulting fragment size.
#'
#'
#' In GBS, large and very small sequences are less likely to be sequenced.
#' This function attempts to replicate that.
#'
#'
#' @param object Either a \code{dna_set} or \code{variants} object that has been
#'     digested (see \code{\link{digest_reference}} and \code{\link{digest_variants}}).
#' @param n_cores Number of cores to use for processing. This argument is ignored if
#'     OpenMP is not enabled. Defaults to 1.
#' @param increment The width of fragment sizes to group together when calculating
#'      the availability probability. Defaults to 1.
#' @param in_place Boolean for whether to edit the object in place.
#'     Defaults to \code{FALSE}.
#'
#' @return If \code{in_place == FALSE}, then a \code{dna_set} or \code{variants} object
#'     that has been filtered. If \code{in_place == TRUE}, \code{NULL}.
#'
#' @export
#'
#' @examples
#'
#' ref_genome <- dna_set$new(rando_seqs(100, mean_len = 1e3, sd_len = 1e2))
#' digest_reference(ref_genome, 'ApeKI', n_cores = 1, in_place = TRUE)
#' size_filter(ref_genome)
#'
#' vars <- make_variants(ref_genome, 2)
#' digest_variants(vars, 'ApeKI', n_cores = 1, in_place = TRUE)
#' size_filter(vars, in_place = TRUE)
#' vars
#'
size_filter <- function(object, n_cores = 1, increment = 1, in_place = FALSE) {

    if (length(object$digests) == 0) stop("Input object should be digested.")

    if (is(object, 'dna_set')) {
        frag_inds <- filter_reference_frags(object$digests, object$sequence_set,
                                            n_cores, increment)
    } else if (is(object, 'variants')) {
        frag_inds <- filter_variants_frags(object$digests, object$reference,
                                           object$variant_set, n_cores, increment)
    } else {
        stop("Input object must be a dna_set or variants object.")
    }
    if (in_place) {
        object$frags_keep <- frag_inds
        invisible(NULL)
    } else {
        out <- object$clone(deep = TRUE)
        out$frags_keep <- frag_inds
        return(out)
    }
}
