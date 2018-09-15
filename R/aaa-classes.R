# ref_genome class----
#' An R6 class representing a reference genome.
#'
#'
#' \emph{Note:} Do NOT change fields in this class directly.
#' It will cause your R session to do bad things.
#' (Ever seen the bomb popup on RStudio? Manually mess with these fields and you
#' surely will.)
#'
#'
#'
#' @field genome An \code{externalptr} to a C++ object storing the sequences
#'     representing the genome.
#' @field digests An \code{externalptr} to a C++ object storing the digestion of the
#'     genome, if a digestion has been carried out. It's \code{NULL} otherwise.
#'
#' @return An object of class \code{ref_genome}.
#'
#' @docType class
#'
#' @seealso \code{\link{read_fasta}} \code{\link{create_genome}}
#'
#' @export
#'
#' @format An \code{\link[R6]{R6Class}} generator object
#'
ref_genome <- R6::R6Class(

    "ref_genome",

    public = list(

        genome = NULL,
        digests = NULL,

        initialize = function(genome_ptr) {
            if (!inherits(genome_ptr, "externalptr")) {
                stop("\nWhen initializing a ref_genome object, you need to use ",
                     "an externalptr object.", call. = FALSE)
            }
            self$genome <- genome_ptr
        },

        print = function(...) {
            print_ref_genome(self$genome)
            invisible(self)
        },

        # -------*
        # Merge all ref_genome genome sequences into one
        # -------*
        merge = function() {
            merge_sequences(self$genome)
            invisible(self)
        },

        # -------*
        # Filter ref_genome sequences by size or for a proportion of total bases
        # -------*
        filter = function(threshold, method = c("size", "prop")) {
            method <- match.arg(method)
            min_seq_size <- 0
            out_seq_prop <- 0
            # Filling in the necessary parameter and checking for sensible inputs
            if (!is.numeric(threshold)) {
                stop("\nWhen filtering ref_genome genome, the threshold must be numeric.",
                     call. = FALSE)
            }
            if (method == "size") {
                if (threshold < 1 | threshold %% 1 != 0) {
                    stop("\nWhen filtering ref_genome genome based on sequence ",
                         "sizes, the threshold must be a whole number greater than 0.",
                         call. = FALSE)
                }
                min_seq_size <- threshold
            } else {
                if (threshold >= 1 | threshold <= 0) {
                    stop("\nWhen filtering ref_genome genome based on a proportion of ",
                         "total bases, the threshold must be > 0 and < 1",
                         call. = FALSE)
                }
                out_seq_prop <- threshold
            }
            filter_sequences(self$genome, min_seq_size, out_seq_prop)
            invisible(self)
        }

    ),

    active = list(

        # -------*
        # Get vector of sequence sizes
        # -------*
        sizes = function() {
            if (!inherits(self$genome, "externalptr")) {
                stop("\nYou're attempting to get sizes using a ref_genome object with ",
                     "a genome field that is not an externalptr. ",
                     "Restart by reading a FASTA file or by simulating a genome. ",
                     "And do NOT change the genome field manually.",
                     call. = FALSE)
            }
            return(view_ref_genome_seq_sizes(self$genome))
        }

    )

)



ref_genome$lock()




# mevo class----
#' An R6 class containing information needed for molecular evolution.
#'
#' You shouldn't need to interact with this class much, if at all.
#' This class is designed to be made in `make_mevo`, then to store information for use
#' in `create_variants`.
#'
#' @field Q A matrix of substitution rates for each nucleotide.
#' @field pi_tcag Vector of nucleotide equilibrium frequencies for "T", "C", "A", and
#'     "G", respectively.
#' @field insertion_rates Vector of insertion rates by length.
#' @field deletion_rates Vector of deletion rates by length.
#' @field gamma_mats List of matrices specifying "gamma distances" (see definition in
#'     `?make_mevo`) for each sequence.
#' @field chunk_size The size of "chunks" of sequences to first sample uniformly
#'     before doing weighted sampling by rates for each sequence location.
#'     See `?make_mevo` for more information.
#'
#'
#' @return An object of class \code{mevo}.
#'
#' @docType class
#'
#' @seealso \code{\link{make_mevo}}
#'
#' @export
#'
#' @format An \code{\link[R6]{R6Class}} generator object
#'
mevo <- R6::R6Class(

    "mevo",

    public = list(

        Q = NULL,
        pi_tcag = NULL,
        insertion_rates = NULL,
        deletion_rates = NULL,
        gamma_mats = NULL,
        chunk_size = NULL,

        initialize = function(sub_info,
                              insertion_rates,
                              deletion_rates,
                              gamma_mats,
                              chunk_size) {

            self$Q <- sub_info$Q
            self$pi_tcag <- sub_info$pi_tcag
            self$insertion_rates <- insertion_rates
            self$deletion_rates <- deletion_rates
            self$gamma_mats <- gamma_mats
            self$chunk_size <- chunk_size

        },


        print = function(digits = max(3, getOption("digits") - 3), ...) {
            fmt <- paste0("%.", digits, "f")
            cat("< Molecular evolution info >\n")

            cat("# Equilibrium densities:\n")
            cat("  ", sprintf(fmt, self$pi_tcag), "\n")

            cat("# Chunk size: ", self$chunk_size, "\n", sep = "")

            cat("# Among-site variability: ")
            using_gammas <- !all(sapply(self$gamma_mats,
                                        function(x) nrow(x) == 1 && all(x[,2] == 1)))
            cat(using_gammas, "\n")

            cat("# Insertion rates:")
            if (length(self$insertion_rates) == 0) {
                cat(" <none>\n")
            } else if (length(self$insertion_rates) > 10) {
                cat("\n")
                cat(sprintf(fmt, self$insertion_rates[1:10]), "...\n")
            } else {
                cat("\n")
                cat(sprintf(fmt, self$insertion_rates), "\n")
            }
            cat("# Deletion rates:")
            if (length(self$deletion_rates) == 0) {
                cat(" <none>\n")
            } else if (length(self$deletion_rates) > 10) {
                cat("\n ", sprintf(fmt, self$deletion_rates[1:10]), "...\n")
            } else {
                cat("\n ", sprintf(fmt, self$deletion_rates), "\n")
            }

            cat("# Substitution rate matrix:\n")
            prmatrix(m$Q, digits = digits,
                     rowlab = paste("  ", c("T", "C", "A", "G")),
                     collab = c("T", "C", "A", "G"))

            invisible(self)
        },


        # -------*
        # Convert to a XPtr<[Chunk]MutationSampler> object
        # -------*
        to_ptr = function() {

            stopifnot(is.numeric(self$chunk_size) & !is.na(self$chunk_size))

            if (self$chunk_size <= 0) {
                sampler_ptr <- make_mutation_sampler_base(self$Q,
                                                          self$pi_tcag,
                                                          self$insertion_rates,
                                                          self$deletion_rates)
            } else {
                sampler_ptr <- make_mutation_sampler_chunk_base(self$Q,
                                                                self$pi_tcag,
                                                                self$insertion_rates,
                                                                self$deletion_rates,
                                                                self$chunk_size)
            }

            return(sampler_ptr)
        },

        mu = function() {
            # Indel rates (same for each nucleotide):
            ins <- sum(self$insertion_rates * 0.25)
            del <- sum(self$deletion_rates * 0.25)
            # Average mutation rate among all nucleotides:
            mu <- sum({rowSums(self$Q) + (ins + del)} * self$pi_tcag)
            return(mu)
        }

    ),


    private = list()

)

mevo$lock()






# variants class----
#' An R6 class representing haploid variants from a reference genome.
#'
#' \emph{Note:} Do NOT change fields in this class directly. It will cause your
#' R session to do bad things.
#' (Ever seen the bomb popup on RStudio? Manually mess with these fields and you
#' surely will.)
#'
#'
#' @field genome An \code{externalptr} to a C++ object storing the sequences
#'     representing the genome.
#' @field digests An \code{externalptr} to a C++ object storing the digestion of the
#'     genome, if a digestion has been carried out. It's \code{NULL} otherwise.
#' @field reference An \code{externalptr} to a C++ object storing the sequences
#'     representing the genome.
#'     There are a few extra notes for this field:
#'     \itemize{
#'         \item \strong{This point is the most important.}
#'             Since it's a pointer, if you make any changes to the reference genome
#'             that it points to, those changes will also show up in the \code{variants}
#'             object. For example, if you make a \code{variants} object \code{V}
#'             based on an existing \code{ref_genome} object \code{R}, then you merge
#'             sequences in \code{R}, \code{V} will now have merged sequences.
#'             If you've already started adding mutations to \code{V},
#'             then all the indexes used to store those mutations will be inaccurate.
#'             So when you do anything with \code{V} later, your R session will crash.
#'         \item This field is private so cannot be accessed directly.
#'         \item If a \code{ref_genome} object is used to create a \code{variants} object,
#'             don't worry about later deleting the \code{ref_genome} object.
#'     }
#'
#'
#' @return An object of class \code{variants}.
#'
#' @docType class
#'
#' @seealso \code{\link{create_variants}}
#'
#' @export
#'
#' @format An \code{\link[R6]{R6Class}} generator object
#'
variants <- R6::R6Class(

    "variants",

    public = list(

        genomes = NULL,
        digests = NULL,

        initialize = function(genomes_ptr, reference_ptr) {
            if (!inherits(genomes_ptr, "externalptr")) {
                stop("\nWhen initializing a variants object, you need to use ",
                     "an externalptr object.", call. = FALSE)
            }
            if (!inherits(reference_ptr, "externalptr")) {
                stop("\nWhen initializing a variants object, you need to use ",
                     "an externalptr object.", call. = FALSE)
            }
            self$genomes <- genomes_ptr
            private$reference <- reference_ptr
        },

        print = function() {
            print_var_set(self$genomes)
            invisible(self)
        }
    ),


    private = list(
        # This should store a `XPtr<RefGenome>` to make sure it doesn't
        # go out of scope:
        reference = NULL
    )

)

variants$lock()

