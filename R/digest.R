#
# Digest genome(s)
#


#' Check validity, adjust for degeneracy, and remove duplicates for input enzymes.
#'
#' @param enzyme_names Name of enzymes that should be present in the `binding_sites`
#'     dataset.
#' @param custom_enzymes Custom enzymes for those not present in internal data.
#'
#' @noRd
#'
process_enzymes <- function(enzyme_names, custom_enzymes) {

    if (missing(enzyme_names) & missing(custom_enzymes)) {
        stop("\nWhen digesting a reference genome, you must provide either an ",
             "enzyme name, a custom enzyme, or both.",
             call. = FALSE)
    }
    enzyme_seqs <- character(0)
    if (!missing(enzyme_names)) {
        enzyme_seqs <- character(length(enzyme_names))
        for (i in 1:length(enzyme_names)) {
            x <- enzyme_names[i]
            z <- binding_sites$sequence[binding_sites == x]
            if (length(z) == 0) {
                stop(paste0("\nEnzyme name \"", x, "\" not found."),
                     call. = TRUE)
            }
            if (grepl("\\(|\\)", z)) {
                stop(paste0("\n", x, " is a non-palindromic enzyme. ",
                            "See the link in ?binding_sites for ",
                            "more info."),
                     call. = TRUE)
            }
            if (!grepl("/", z)) {
                stop(paste("\n", x, "doesn't have a cleavage site",
                           "according to NEB.",
                           "Choose another enzyme,",
                           "or provide a custom enzyme if you're",
                           "sure I've erred in this error."),
                     call. = TRUE)
            }
            enzyme_seqs[i] <- z
        }

    }
    if (!missing(custom_enzymes)) {
        stopifnot(inherits(custom_enzymes, "character"))

        allowed_chars <- sort(unique(c(nucleobase_legend, "/", recursive = TRUE)))
        allowed_chars <- paste(allowed_chars, collapse = "")
        weird_chars <- grepl(sprintf("[^%s]", allowed_chars), custom_enzymes)
        if (any(weird_chars)) {
            stop(paste("\nThe only allowed characters for restriction enzymes are",
                       allowed_chars),
                 call. = FALSE)
        }

        no_cleavage <- ! grepl("/", custom_enzymes)

        if (any(no_cleavage)) {
            stop("\nYou must include a cleavage site for each custom restriction ",
                 "site using \"/\" (e.g., \"TTA/TAA\" for the enzyme AanI).",
                 call. = TRUE)
        }

        enzyme_seqs <- c(enzyme_seqs, custom_enzymes)
    }

    # Get number of nucleotides before the cleavage site:
    len5s <- get_precleavage_lens(enzyme_seqs)

    # Now we no longer want/need the cleavage indicators:
    bind_sites <- gsub("/", "", enzyme_seqs)

    bind_sites <- expand_seqs(bind_sites)

    return(list(len5s = len5s, bind_sites = bind_sites))
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
#' \dontrun{
#'
#' ref_genome <- dna_set$new(gemino:::rando_seqs(100, mean_len = 1e3, sd_len = 1e2))
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
#' }
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

    if (inherits(object, 'dna_set')) {
        digest_reference(object, enzyme_names, n_cores, enz_list, in_place)
    } else if (inherits(object, 'variants')) {
        digest_variants(object, enzyme_names, n_cores, chunk_size, enz_list, in_place)
    } else {
        stop("Input object must be a dna_set or variants object.")
    }
    invisible(NULL)
}
