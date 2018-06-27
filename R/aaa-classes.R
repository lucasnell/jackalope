#' An R6 class representing a reference genome
#'
#' Note: all fields are private so that they cannot be manipulated manually.
#' Trust me that this is a good thing.
#'
#' @field genome An \code{externalptr} to a C++ object storing the sequences
#'     representing the genome.
#' @field digests An \code{externalptr} to a C++ object storing the digestion of the
#'     genome, if a digestion has been carried out. It's \code{NULL} otherwise.
#'
#' @return An object of class \code{reference}.
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @export
#'
#' @format An \code{\link[R6]{R6Class}} generator object
#'
reference <- R6::R6Class(

    "reference",

    public = list(

        initialize = function(genome_ptr) {
            if (!inherits(genome_ptr, "externalptr")) {
                stop("\nWhen initializing a reference object, you need to use ",
                     "an externalptr object.", call. = FALSE)
            }
            private$genome <- genome_ptr
        },

        print = function(...) {
            "print reference object"
            print_ref_genome(genome)
        },

        digest = function(bind_sites,
                          len5s,
                          chunk_size = 0,
                          n_cores = 1) {
            "Digest reference genome based on restriction enzyme(s)"
            # Warning should be removed and "digest_ref(...)" lines should be
            # uncommented when digest returns a C++ class:
            warning("Digestion not yet implemented for reference class", call. = FALSE)
            # private$digests <- digest_ref(genome, bind_sites, len5s, chunk_size, n_cores)
            invisible(self)
        },

        merge = function() {
            "Merge all reference genome sequences into one"
            merge_sequences(genome)
            invisible(self)
        },

        filter = function(threshold, method = c("size", "prop")) {
            "Filter reference genome sequences by size or for a proportion of total bases"
            method <- match.arg(method)
            min_seq_size <- 0
            out_seq_prop <- 0
            # Filling in the necessary parameter and checking for sensible inputs
            if (!is.numeric(threshold)) {
                stop("\nWhen filtering reference genome, the threshold must be numeric.",
                     call. = FALSE)
            }
            if (method == "size") {
                if (threshold < 1 | threshold %% 1 != 0) {
                    stop("\nWhen filtering reference genome based on sequence ",
                         "sizes, the threshold must be a whole number greater than 0.",
                         call. = FALSE)
                }
                min_seq_size <- threshold
            } else {
                if (threshold >= 1 | threshold <= 0) {
                    stop("\nWhen filtering reference genome based on a proportion of ",
                         "total bases, the threshold must be > 0 and < 1",
                         call. = FALSE)
                }
                out_seq_prop <- threshold
            }
            filter_sequences(genome, min_seq_size, out_seq_prop)
            invisible(self)
        }

    ),

    private = list(
        genome = NULL,
        digests = NULL
    )
)




#' An R6 class representing haploid variants from a reference genome
#'
#' Note: all fields are private so that they cannot be manipulated manually.
#' Trust me that this is a good thing.
#'
#'
#' @field genome An \code{externalptr} to a C++ object storing the sequences
#'     representing the genome.
#' @field digests An \code{externalptr} to a C++ object storing the digestion of the
#'     genome, if a digestion has been carried out. It's \code{NULL} otherwise.
#'
#' @return An object of class \code{variants}.
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @export
#'
#' @format An \code{\link[R6]{R6Class}} generator object
#'
variants <- R6::R6Class(

    "variants",

    public = list(

        initialize = function(genomes_ptr) {
            if (!inherits(genomes_ptr, "externalptr")) {
                stop("\nWhen initializing a variants object, you need to use ",
                     "an externalptr object.", call. = FALSE)
            }
            private$genomes <- genomes_ptr
        },

        print = function() {
            "print reference object"
            print_var_set(genomes)
        },

        digest = function(bind_sites,
                          len5s,
                          chunk_size = 1000,
                          n_cores = 1) {
            "Digest reference genome based on restriction enzyme(s)"
            # Warning should be removed and "digest_ref(...)" lines should be
            # uncommented when digest returns a C++ class:
            warning("Digestion not yet implemented for variants class", call. = FALSE)
            # private$digests <- digest_var_set(genomes, bind_sites, len5s,
            #                                   chunk_size, n_cores)
            invisible(self)
        }
    ),

    private = list(
        genomes = NULL,
        digests = NULL
    )
)

