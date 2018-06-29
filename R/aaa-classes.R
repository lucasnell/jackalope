#' An R6 class representing a reference genome
#'
#' \emph{Note:} Do NOT change fields in this class directly. It will cause your
#' R session to do bad things.
#' (Ever seen the bomb popup on RStudio? Manually mess with these fields and you
#' surely will.)
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

        genome = NULL,
        digests = NULL,

        initialize = function(genome_ptr) {
            if (!inherits(genome_ptr, "externalptr")) {
                stop("\nWhen initializing a reference object, you need to use ",
                     "an externalptr object.", call. = FALSE)
            }
            self$genome <- genome_ptr
        },

        print = function(...) {
            "print reference object"
            print_ref_genome(self$genome)
            invisible(self)
        },

        merge = function() {
            "Merge all reference genome sequences into one"
            merge_sequences(self$genome)
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
            filter_sequences(self$genome, min_seq_size, out_seq_prop)
            invisible(self)
        }

    )

)



reference$lock()




#' An R6 class representing haploid variants from a reference genome
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
#'             based on an existing \code{reference} object \code{R}, then you merge
#'             sequences in \code{R}, \code{V} will now have merged sequences.
#'             If you've already started adding mutations to \code{V},
#'             then all the indexes used to store those mutations will be inaccurate.
#'             So when you do anything with \code{V} later, your R session will crash.
#'         \item This field is private so cannot be accessed directly.
#'         \item If a \code{reference} object is used to create a \code{variants} object,
#'             don't worry about later deleting the \code{reference} object.
#'     }
#'
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

        genomes = NULL,
        digests = NULL,

        initialize = function(genomes_ptr, reference_ptr) {
            if (!inherits(genomes_ptr, "externalptr")) {
                stop("\nWhen initializing a variants object, you need to use ",
                     "an externalptr object.", call. = FALSE)
            }
            self$genomes <- genomes_ptr
            private$reference <- reference_ptr
        },

        print = function() {
            "print reference object"
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

