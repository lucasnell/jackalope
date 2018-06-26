
#'
#' An R6 class to represent a set of DNA sequences
#'
#'
#'
#' @field sequences A character vector of nucleotide sequences.
#' @field seq_names A character vector of sequence names. Defaults to
#'     \code{paste0('seq_', 1:length(sequences))}.
#'
#' @return An object of class dna_set.
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dna_set$new(c("AACT", "TGGCA", "AAAAAATTTTTT"), LETTERS[1:3])
#' }
#'
#' @format An \code{\link[R6]{R6Class}} generator object
#'
dna_set <- R6::R6Class(
    "dna_set",
    public = list(

        sequence_set = NULL,
        digests = list(),
        frags_keep = list(),
        binding_sites = NULL,

        initialize = function(sequences = NULL, seq_names = NULL,
                              fasta_file = NULL, fai_file = NULL,
                              cut_names = TRUE, rm_soft_mask = TRUE) {

            if (!is.null(fasta_file) & is.null(sequences)) {
                stopifnot(is.character(fasta_file), length(fasta_file) == 1,
                          is.logical(cut_names), is.logical(rm_soft_mask),
                          length(cut_names) == 1, length(rm_soft_mask) == 1)
                if (grepl('.fa$|.fasta$|fa.gz$|fasta.gz$', fasta_file)) {
                    self$sequence_set <- read_fasta(fasta_file, fai_file,
                                                    cut_names, rm_soft_mask)
                } else {
                    stop(paste("The fasta filename is ambiguous. Please use either",
                               "the .fa or .fasta extension."))
                }
            } else if (!is.null(sequences) & is.null(fasta_file)) {
                if (is.null(seq_names)) {
                    seq_names <- paste0('seq_', 1:length(sequences))
                }
                stopifnot(is.character(sequences), is.character(seq_names),
                          length(sequences) == length(seq_names))
                self$sequence_set <- SequenceSet_characters(sequences, seq_names)
            } else {
                stop(paste("Either specify a fasta file OR sequences. Not both."))
            }
        },

        print = function(...) {
            SummarizeSequenceSet(self$sequence_set, options('width')$width)
        },

        merge = function(n_merged_scaffs) {
            merge_scaffolds(self$sequence_set, n_merged_scaffs)
        },

        names = function() {
            seq_names_SequenceSet(self$sequence_set)
        },

        gc_content = function(seq_name = NULL, start = NULL, stop = NULL) {
            if (is.null(seq_name)) {  # get all scaffolds' GC content
                gc <- all_scaff_gc_SequenceSet(self$sequence_set)
            } else if (is.null(start) | is.null(stop)) {  # get one scaffold's GC content
                gc <- one_scaff_gc_SequenceSet(self$sequence_set, seq_name)
            } else if (length(start) == 1) {  # get one range's GC content
                gc <- one_range_gc_SequenceSet(self$sequence_set, seq_name, start, stop)
            } else {  # get multiple ranges' GC content
                gc <- mult_ranges_gc_SequenceSet(self$sequence_set, seq_name, start, stop)
            }
            return(gc)
        },

        filter = function(min_size = NULL, prop = NULL) {
            if ((is.null(min_size) & is.null(prop)) |
                (!is.null(min_size) & !is.null(prop))) {
                stop(paste("\nProvide exactly one of the following to filter:\n",
                           "(1) minimum scaffold size (in bp)\n",
                           "(2) minimum proportion of total genomic sequence to retain"))
            }
            if (is.null(min_size)) min_size <- 0
            if (is.null(prop)) prop <- 0
            stopifnot(prop < 1, prop >= 0, min_size >= 1 | min_size == 0)
            filter_scaffolds(self$sequence_set, out_scaff_prop = prop,
                             out_scaff_size = min_size)
        }

    ),

    private = list(
        # With x$clone(deep=TRUE) is called, the deep_clone gets invoked once for
        # each field, with the name and value.
        deep_clone = function(name, value) {
            if (name == "sequence_set") {
                # Make a copy of this C++ object
                SequenceSet_copy(value)
            } else {
                value
            }
        }
    )
)

















#' Genome variants class.
#'
#' An R6 class to represent variants from a reference genome.
#'
#'
#' @field reference An external pointer (to a C++ SequenceSet) representing the
#'     reference genome. An external pointer (to a C++ SequenceSet),
#'     \code{\link{dna_set}}, or character vector can be used for initialization.
#' @field nucleos A list of character vectors, one vector for each variant.
#'     Each string in a vector represents the nucleotides present at the segregating
#'     sites for a single scaffold from a single variant.
#'     See the examples below for more.
#' @field sites A list of integer vectors, one vector for each scaffold.
#'     Each vector contains the locations (1-indexed) for each segregating site on
#'     that scaffold.
#' @field digests A list of lists, each sub-list containing multiple vectors representing
#'     the locations of cut sites for a given variant on a given scaffold.
#'     Indexing this list would be done as such:
#'     \code{var_obj$digests[[variant_index]][[scaffold_index]][position_index]}.
#'     Cut site positions represent starting positions but are indexed using C++'s
#'     0-based indexing.
#'     For example, if the position column's first value is 10,
#'     this means that the first fragment from that scaffold ends at the 10th position.
#'     The 10th position is indexed using \code{[10]} in R, but using \code{[9]} in C++.
#'     For each scaffold, there will be \eqn{n + 1} resulting fragments,
#'     where \eqn{n} is the number of items in the vector for that scaffold.
#'
#' @return An object of class \code{variants}.
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @usage
#' variants$new(reference, nucleos, sites, digest_common = list(),
#'              digest_unique = list(), binding_sites = NULL)
#'
#' @export
#'
#' @format An \code{\link[R6]{R6Class}} generator object
#'
variants <- R6::R6Class(
    "variants",
    public = list(
        reference = NULL,
        variant_set = NULL,
        digests = list(),
        frags_keep = list(),
        binding_sites = NULL,
        initialize = function(reference, variant_set,
                              digests = list(),
                              binding_sites = NULL){
            if (is(reference, 'dna_set')) {
                self$reference <- reference$sequence_set
            } else if (is(reference, 'externalptr')) {
                self$reference <- reference
            } else if (is(reference, 'character')) {
                self$reference <- SequenceSet_characters(
                    reference, paste0('seq_', 1:length(reference)))
            } else {
                stop(paste("Reference must be an external pointer,",
                           "dna_set, or character."))
            }
            if (is(variant_set, 'externalptr')) {
                self$variant_set <- variant_set
            } else stop("variant_set must be an external pointer.")

            self$digests <- digests
            self$binding_sites <- binding_sites
        },
        get_variant = function(variant_num) {
            var_scaffs <- variants_retr_var(variant_num, self$reference,
                                       self$variant_set)
            return(dna_set$new(var_scaffs))
        },
        get_scaff = function(scaff_num, variant_num) {
            scaff <- variants_retr_scaff(scaff_num, variant_num,
                                    self$reference, self$variant_set)
            return(scaff)
        },
        get_seq = function(start, stop, scaff_num, variant_num) {
            length_out <- stop - start + 1
            seq <- variants_retr_seq(start, length_out,
                                scaff_num, variant_num,
                                self$reference, self$variant_set)
            return(seq)
        },
        print = function(...) {
            # Variants info
            cat(cpp_merge_str(rep(' ', ceiling({options('width')$width - 21}/2))),
                '<< Variants object >>\n')

            SummarizeVariantSet(self$variant_set)

            # Digestion info
            digested <- !is.null(self$binding_sites)
            cat('\n* Digested:', digested, '\n')
            if (digested) {
                prime5_prime3s <- sapply(
                    seq(1, length(self$binding_sites), 2),
                    function(i) cpp_merge_str(self$binding_sites[i:(i+1)]))
                cat('* Binding sites: ')
                cat(paste0(prime5_prime3s, collapse = ' '), '\n')
            }

            # Reference genome info
            cat('\n', cpp_merge_str(rep(' ', ceiling({options('width')$width - 28}/2))),
                '<< Reference genome info: >>\n')
            SummarizeSequenceSet(self$reference, options('width')$width)
        }
    ),
    # These active bindings allow access to these fields quickly and the same as if they
    # were a field in this class (even though they're present in the C++ class VariantSet)
    active = list(
        n_variants = function() {
            return(n_variants_VS(self$variant_set))
        },
        total_segr_sites = function() {
            return(total_segr_sites_VS(self$variant_set))
        },
        nucleos = function(variant, scaff) {
            return(nucleos_VS(self$variant_set, variant - 1, scaff - 1))
        },
        sites = function(variant, scaff) {
            return(sites_VS(self$variant_set, variant - 1, scaff - 1))
        },
        scaffold_lengths = function(variant, scaff) {
            return(scaffold_length_VS(self$variant_set, variant - 1, scaff - 1))
        }
    ),

    private = list(
        # With x$clone(deep=TRUE) is called, the deep_clone gets invoked once for
        # each field, with the name and value.
        deep_clone = function(name, value) {
            if (name == "variant_set") {
                # Make a copy of this C++ object
                VariantSet_copy(value)
            } else if (name == "reference") {
                # Make a copy of this C++ object
                SequenceSet_copy(value)
            } else {
                value
            }
        }
    )
)

