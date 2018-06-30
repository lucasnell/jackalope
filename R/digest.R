#
# Digest genome(s)
#



#' Check validity, adjust for degeneracy, and remove duplicates for input enzymes.
#'
#' @inheritParams digest
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
            z <- gemino::binding_sites$sequence[gemino::binding_sites == x]
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

        allowed_chars <- sort(unique(c(gemino::nucleobase_legend, "/",
                                       recursive = TRUE)))
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







#' Digest reference genome or variants from a genome.
#'
#' \emph{Note:} This will override any digestions currently in place in the
#' object. If you want to add a new digestion, re-run this function with the names
#' of all enzymes you're interested in included in the \code{enzyme_names} argument.
#'
#'
#' @inheritParams create_genome
#' @param object Either a \code{ref_genome} or \code{variants} object.
#' @param enzyme_names Name of enzymes that should be present inside the package data
#'     object `binding_sites` (in the `enzyme` column).
#' @param custom_enzymes Custom enzymes for those not present in internal data.
#' @param chunk_size The size of chunks to break up sqeuences into when digesting
#'     a \code{ref_genome} or \code{variants} object.
#'     Changing this might affect performance, for better or worse.
#'     The default worked best on my computer. Defaults to \code{1000}.
#'
#' @return All changes occur in place, but the input object is returned invisibly
#'     so that piping works.
#'
#'
#' @export
#'
#'
digest <- function(object,
                   enzyme_names,
                   custom_enzymes,
                   chunk_size = 0,
                   n_cores = 1) {

    if (missing(enzyme_names) & missing(custom_enzymes)) {
        stop("\nWhen digesting a reference genome or variants, you must provide ",
             "either an enzyme name, a custom enzyme, or both.",
             call. = FALSE)
    }

    enz_info <- process_enzymes(enzyme_names, custom_enzymes)

    if (inherits(object, 'ref_genome')) {
        if (!inherits(object$genome, "externalptr")) {
            stop("\nYou're attempting a digestion on a ref_genome object with ",
                 "a genome field that is not an externalptr. ",
                 "Restart by reading a FASTA file or by simulating a genome. ",
                 "And do NOT change the genome field manually.",
                 call. = FALSE)
        }
        object$digests <- digest_ref(object$genome,
                                     enz_info$bind_sites, enz_info$len5s,
                                     chunk_size, n_cores)
    } else if (inherits(object, 'variants')) {
        if (!inherits(object$genomes, "externalptr")) {
            stop("\nYou're attempting a digestion on a variants object with ",
                 "a genomes field that is not an externalptr. ",
                 "Restart by reading a FASTA file or by simulating a genome, ",
                 "then generating variants. ",
                 "And do NOT change the genomes field manually.",
                 call. = FALSE)
        }
        object$digests <- digest_var_set(object$genomes,
                                         enz_info$bind_sites, enz_info$len5s,
                                         chunk_size, n_cores)
    } else {
        stop("\nA digestion can only be applied to a ref_genome or variants object.",
             call. = FALSE)
    }
    invisible(object)
}

