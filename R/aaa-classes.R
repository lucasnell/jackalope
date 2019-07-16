# >> CLASS ref_genome----
#' An R6 class representing a reference genome.
#'
#'
#' \emph{Note:} This class wraps a pointer to a C++ object, so
#' do NOT change fields in this class directly.
#' It will cause your R session to do bad things.
#' (Ever seen the bomb popup on RStudio? Manually mess with these fields and you
#' surely will.)
#' For safe ways of manipulating the reference genome, see the "Methods" section.
#'
#'
#' @field genome An \code{externalptr} to a C++ object storing the sequences
#'     representing the genome.
#'
#' @section Methods:
#' \strong{Viewing information:}
#' \describe{
#'     \item{`n_seqs()`}{View the number of sequences.}
#'     \item{`sizes()`}{View vector of sequence sizes.}
#'     \item{`names()`}{View vector of sequence names.}
#'     \item{`sequence(seq_ind)`}{View a sequence string based on an index,
#'         `seq_ind`.}
#'     \item{`gc_prop(seq_ind, start, end)`}{View the GC proportion for a range within a
#'         reference sequence.}
#'     \item{`nt_prop(nt, seq_ind, start, end)`}{View the proportion of a range within a
#'         reference sequence that is of nucleotide `nt`.}
#' }
#' \strong{Editing information:}
#' \describe{
#'     \item{`set_names(new_names)`}{Set names for all sequences.
#'         `new_names` is a character vector of what to change names to, and it must
#'         be the same length as the # sequences.}
#'     \item{`clean_names()`}{Clean sequence names, converting `" :;=%,\\|/\"\'"`
#'         to `"_"`.}
#'     \item{`add_seqs(new_seqs, new_names = NULL)`}{Add one or more sequences
#'         directly. They can optionally be named (using `new_names`).
#'         Otherwise, their names are auto-generated.}
#'     \item{`rm_seqs(seq_names)`}{Remove one or more sequences based on names in
#'         the `seq_names` vector.}
#'     \item{`merge_seqs()`}{Merge all sequences into one after first shuffling
#'         their order.}
#'     \item{`filter_seqs(threshold, method)`}{Filter sequences by size
#'         (`method = "size"`) or for a proportion of total bases (`method = "prop"`).
#'         For the latter, sequences are first size-sorted, then the largest `N`
#'         sequences are retained that allow at least
#'         `threshold * sum(<all sequence sizes>)` base pairs remaining after
#'         filtering.}
#'     \item{`replace_Ns(pi_tcag, n_threads = 1, show_progress = FALSE)`}{Replace
#'         `N`s in reference sequence with nucleotides sampled with probabilities
#'         given in `pi_tcag`.
#'         You can optionally use multiple threads (`n_threads` argument) and/or
#'         show a progress bar (`show_progress`).}
#' }
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
#' @importFrom R6 R6Class
#'
ref_genome <- R6Class(

    "ref_genome",

    public = list(

        genome = NULL,

        initialize = function(genome_ptr) {
            if (!inherits(genome_ptr, "externalptr")) {
                stop("\nWhen initializing a ref_genome object, you need to use ",
                     "an externalptr object.", call. = FALSE)
            }
            self$genome <- genome_ptr
        },

        print = function(...) {
            private$check_ptr()
            print_ref_genome(self$genome)
            invisible(self)
        },

        # ----------*
        # __view__ ----
        # ----------*
        # Get # sequences
        n_seqs = function() {
            private$check_ptr()
            return(view_ref_genome_nseqs(self$genome))
        },

        # Get vector of sequence sizes
        sizes = function() {
            private$check_ptr()
            return(view_ref_genome_seq_sizes(self$genome))
        },

        # Get vector of sequence names
        names = function() {
            private$check_ptr()
            return(view_ref_genome_seq_names(self$genome))
        },

        # Extract one reference sequence
        sequence = function(seq_ind) {
            private$check_ptr()
            if (!single_integer(seq_ind, 1, self$n_seqs())) {
                err_msg("sequence", "seq_ind", "integer in range [1, <# sequences>]")
            }
            return(view_ref_genome_seq(self$genome, seq_ind - 1))
        },
        # GC proportion for part of one reference sequence
        gc_prop = function(seq_ind, start, end) {
            private$check_pos(seq_ind, start, "nt_prop", "start")
            private$check_pos(seq_ind, end, "nt_prop", "end")
            if (end < start) err_msg("gc_prop", "end", ">= `start` arg")
            gcp <- view_ref_genome_gc_content(self$genome, seq_ind - 1,
                                              start - 1, end - 1)
            return(gcp)
        },
        # Nucleotide content for part of one reference sequence
        nt_prop = function(nt, seq_ind, start, end) {
            private$check_pos(seq_ind, start, "nt_prop", "start")
            private$check_pos(seq_ind, end, "nt_prop", "end")
            if (end < start)err_msg("nt_prop", "end", ">= `start` arg")
            if (!is_type(nt, "character", 1) || nchar(nt) != 1) {
                err_msg("nt_prop", "nt", "a single character")
            }
            ntp <- view_ref_genome_nt_content(self$genome, nt, seq_ind - 1,
                                              start - 1, end - 1)
            return(ntp)
        },

        # ----------*
        # __edit__ ----
        # ----------*
        # Change sequence names
        set_names = function(new_names) {
            private$check_ptr()
            if (!is_type(new_names, "character", self$n_seqs())) {
                err_msg("set_names", "new_names", "the same length as # sequences")
            }
            seq_inds <- 0:(length(new_names) - 1)
            set_ref_genome_seq_names(self$genome, seq_inds, new_names)
            invisible(self)
        },

        # Clean sequence names, converting " :;=%,\\|/\"\'" to "_"
        clean_names = function() {
            private$check_ptr()
            clean_ref_genome_seq_names(self$genome)
            invisible(self)
        },

        # Add one or more sequences
        add_seqs = function(new_seqs, new_names = NULL) {
            private$check_ptr()
            self_names <- self$names()
            if (!is_type(new_seqs, "character")) {
                err_msg("rm_seqs","new_seqs", "a character vector")
            }
            if (!is.null(new_names) && !is_type(new_names, "character")) {
                err_msg("rm_seqs","new_names", "NULL or a character vector")
            }
            if (is.null(new_names)) {
                N <- suppressWarnings(as.integer(gsub("seq", "", self$names())))
                N <- if (any(!is.na(N))) 1 + max(N, na.rm = TRUE) else self$n_seqs()
                N <- N:(N+length(new_seqs)-1)
                new_names <- sprintf("seq%i", N)
            }

            if (length(new_names) != length(new_seqs)) {
                err_msg("rm_seqs","new_names", "the same length as `new_seqs`")
            }
            if (anyDuplicated(new_names) != 0) {
                err_msg("rm_seqs","new_names", "a vector with no duplicates")
            }
            if (any(new_names %in% self$names())) {
                err_msg("rm_seqs", "new_names", "a vector containing no elements",
                        "already present as a sequence name in the reference genome")
            }

            add_ref_genome_seqs(self$genome, new_seqs, new_names)
            invisible(self)
        },

        # Remove one or more sequences by name
        rm_seqs = function(seq_names) {
            private$check_ptr()
            self_names <- self$names()
            if (!is_type(seq_names, "character")) {
                err_msg("rm_seqs", "seq_names", "a character vector")
            }
            if (!all(seq_names %in% self_names)) {
                err_msg("rm_seqs", "seq_names", "a character vector of sequence names",
                        "present in the `ref_genome` object. One or more of the names",
                        "provided weren't found")
            }
            if (anyDuplicated(seq_names) != 0) {
                err_msg("rm_seqs", "seq_names", "a vector of *non-duplicated* names")
            }
            seq_inds <- match(seq_names, self_names) - 1
            remove_ref_genome_seqs(self$genome, seq_inds)
            invisible(self)
        },

        # Merge all ref_genome genome sequences into one
        merge_seqs = function() {
            private$check_ptr()
            merge_sequences_cpp(self$genome)
            invisible(self)
        },

        # Filter ref_genome sequences by size or for a proportion of total bases
        filter_seqs = function(threshold, method) {
            private$check_ptr()
            method <- match.arg(method, c("size", "prop"))
            min_seq_size <- 0
            out_seq_prop <- 0
            # Filling in the necessary parameter and checking for sensible inputs
            if (!single_number(threshold)) {
                err_msg("filter_seqs", "threshold", "a single number")
            }
            if (method == "size") {
                if (!single_integer(threshold, 1)) {
                    err_msg("filter_seqs", "threshold", "a single integer >= 1 if",
                            "filtering based on sequence sizes")
                }
                min_seq_size <- threshold
            } else {
                if (!single_number(threshold) || threshold >= 1 || threshold <= 0) {
                    err_msg("filter_seqs", "threshold", "a single number > 0 and < 1 if",
                            "filtering based on a proportion of total bases")
                }
                out_seq_prop <- threshold
            }
            filter_sequences_cpp(self$genome, min_seq_size, out_seq_prop)
            invisible(self)
        },


        replace_Ns = function(pi_tcag,
                              n_threads = 1,
                              show_progress = FALSE) {
            private$check_ptr()
            if (!is_type(pi_tcag, "numeric", 4) || any(pi_tcag < 0) ||
                all(pi_tcag == 0)) {
                err_msg("replace_Ns", "pi_tcag", "a numeric vector of length 4, where",
                        "no values can be < 0 and at least one value must be > 0")
            }
            if (!single_integer(n_threads, 1)) {
                err_msg("replace_Ns", "n_threads", "a single integer >= 1")
            }
            if (!is_type(show_progress, "logical", 1)) {
                err_msg("replace_Ns", "show_progress", "a single logical")
            }

            replace_Ns_cpp(self$genome, pi_tcag, n_threads, show_progress)

            invisible(self)

        }


    ),

    # __private__ -----
    private = list(

        check_ptr = function() {
            if (!inherits(self$genome, "externalptr")) {
                stop("\nSomehow the `genome` field of this `ref_genome` object ",
                     "has been converted to something other than an `externalptr` ",
                     "object. Re-make this `ref_genome` object, don't edit the ",
                     "`genome` field directly, and try again.", call. = FALSE)
            }
        },

        check_pos = function(seq_ind, pos, .fun, .pos_name) {
            private$check_ptr()
            if (!single_integer(seq_ind, 1, self$n_seqs())) {
                err_msg(.fun, "seq_ind", "integer in range [1, <# sequences>]")
            }
            if (!single_integer(pos, 1, self$sizes()[seq_ind])) {
                err_msg(.fun, .pos_name, "integer in range [1, <sequence size>]")
            }
        }

    )

)



ref_genome$lock()







# >> CLASS variants----
#' An R6 class representing haploid variants from a reference genome.
#'
#' \emph{Note:} This class wraps a pointer to a C++ object, so
#' do NOT change fields in this class directly.
#' It will cause your R session to do bad things.
#' (Ever seen the bomb popup on RStudio? Manually mess with these fields and you
#' surely will.)
#' For safe ways of manipulating the variants' information, see the "Methods" section.
#'
#' @field genome An \code{externalptr} to a C++ object storing the sequences
#'     representing the genome.
#' @field reference An \code{externalptr} to a C++ object storing the sequences
#'     representing the genome.
#'     This field is private, so you can't view it, but I'm listing it here
#'     so that I can provide a few extra notes about it:
#'     \itemize{
#'         \item \strong{This point is the most important.}
#'             Since it's a pointer, if you make any changes to the reference genome
#'             that it points to, those changes will also show up in the \code{variants}
#'             object. For example, if you make a \code{variants} object named \code{V}
#'             based on an existing \code{ref_genome} object named \code{R},
#'             then you merge sequences in \code{R},
#'             \code{V} will now have merged sequences.
#'             If you've already started adding mutations to \code{V},
#'             then all the indexes used to store those mutations will be inaccurate.
#'             So when you do anything with \code{V} later, your R session will crash
#'             or have errors.
#'         \item If a \code{ref_genome} object is used to create a \code{variants}
#'             object, deleting the \code{ref_genome} object won't cause issues with
#'             the \code{variants} object.
#'             However, the \code{variants} class doesn't provide methods to edit
#'             sequences, so only remove the \code{ref_genome} object when you're done
#'             editing the reference genome.
#'     }
#'
#'
#' @section Methods:
#' \strong{Viewing information:}
#' \describe{
#'     \item{`n_seqs()`}{View the number of sequences.}
#'     \item{`n_vars()`}{View the number of variants.}
#'     \item{`sizes(var_ind)`}{View vector of sequence sizes for a given variant.}
#'     \item{`seq_names()`}{View vector of sequence names.}
#'     \item{`var_names()`}{View vector of variant names.}
#'     \item{`sequence(var_ind, seq_ind)`}{View a sequence string based on
#'         indices for the sequence (`seq_ind`) and variant (`var_ind`).}
#'     \item{`gc_prop(var_ind, seq_ind, start, end)`}{View the GC proportion for a range
#'         within a variant sequence.}
#'     \item{`nt_prop(nt, var_ind, seq_ind, start, end)`}{View the proportion of a range
#'         within a variant sequence that is of nucleotide `nt`.}
#' }
#' \strong{Editing information:}
#' \describe{
#'     \item{`set_names(new_names)`}{Set names for all variants.
#'         `new_names` is a character vector of what to change names to, and it must
#'         be the same length as the # variants.}
#'     \item{`add_vars(new_names)`}{Add new, named variant(s) to the object.
#'         These variants will have no mutations. If you want to add new variants with
#'         mutations, either re-run `create_variants` or use the `dup_vars` method to
#'         duplicate existing variants.}
#'     \item{`dup_vars(var_names, new_names = NULL)`}{Duplicate existing variant(s) based
#'         on their name(s). You can optionally specify the names of the duplicates
#'         (using `new_names`).
#'         Otherwise, their names are auto-generated.}
#'     \item{`rm_vars(var_names)`}{Remove one or more variants based on names in
#'         the `var_names` vector.}
#'     \item{`add_sub(var_ind, seq_ind, pos, nt)`}{Manually add a substitution
#'         for a given variant (`var_ind`), sequence (`seq_ind`), and position (`pos`).
#'         The reference nucleotide will be changed to `nt`, which should be a single
#'         character.}
#'     \item{`add_ins(var_ind, seq_ind, pos, nts)`}{Manually add an insertion
#'         for a given variant (`var_ind`), sequence (`seq_ind`), and position (`pos`).
#'         The nucleotide(s) `nts` will be inserted after the designated position.}
#'     \item{`add_del(var_ind, seq_ind, pos, n_nts)`}{Manually add a deletion
#'         for a given variant (`var_ind`), sequence (`seq_ind`), and position (`pos`).
#'         The designated number of nucleotides to delete (`n_nts`) will be deleted
#'         starting at `pos`, unless `pos` is near the sequence end and doesn't have
#'         `n_nts` nucleotides to remove; it simply stops at the sequence end in
#'         this case.}
#' }
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
#' @importFrom R6 R6Class
#'
variants <- R6Class(

    "variants",

    public = list(

        genomes = NULL,

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
            private$check_ptr()
            print_var_set(self$genomes)
            invisible(self)
        },


        # ----------*
        # __view__ ----
        # ----------*
        # Get # sequences
        n_seqs = function() {
            private$check_ptr()
            return(view_var_set_nseqs(self$genomes))
        },

        # Get # variants
        n_vars = function() {
            private$check_ptr()
            return(view_var_set_nvars(self$genomes))
        },

        # Get vector of sequence sizes for one variant
        sizes = function(var_ind) {
            private$check_ptr()
            private$check_var_ind(var_ind, "sizes")
            return(view_var_genome_seq_sizes(self$genomes, var_ind - 1))
        },

        # Get vector of sequence names
        seq_names = function() {
            stopifnot(inherits(private$reference, "externalptr"))
            return(view_ref_genome_seq_names(private$reference))
        },

        # Get vector of variant names
        var_names = function() {
            private$check_ptr()
            return(view_var_set_var_names(self$genomes))
        },

        # Extract one variant sequence
        sequence = function(var_ind, seq_ind) {
            private$check_ptr()
            private$check_var_ind(var_ind, "sequence")
            private$check_seq_ind(seq_ind, "sequence")
            return(view_var_genome_seq(self$genomes, var_ind - 1, seq_ind - 1))
        },

        # GC proportion for part of one variant sequence
        gc_prop = function(var_ind, seq_ind, start, end) {
            private$check_pos(var_ind, seq_ind, start, "gc_prop", "start")
            private$check_pos(var_ind, seq_ind, end, "gc_prop", "end")
            if (end < start) err_msg("gc_prop", "end", ">= `start` arg")
            gcp <- view_var_set_gc_content(self$genomes, seq_ind - 1, var_ind - 1,
                                           start - 1, end - 1)
            return(gcp)
        },

        # Nucleotide content for part of one reference sequence
        nt_prop = function(nt, var_ind, seq_ind, start, end) {
            private$check_pos(var_ind, seq_ind, start, "nt_prop", "start")
            private$check_pos(var_ind, seq_ind, end, "nt_prop", "end")
            if (end < start) err_msg("nt_prop", "end", ">= `start` arg")
            if (!is_type(nt, "character", 1) || nchar(nt) != 1) {
                err_msg("nt_prop", "nt", "a single character")
            }
            ntp <- view_var_set_nt_content(self$genomes, nt, seq_ind - 1, var_ind - 1,
                                           start - 1, end - 1)
            return(ntp)
        },


        # ----------*
        # __edit__ ----
        # ----------*

        # Change variant names
        set_names = function(new_names) {
            private$check_ptr()
            if (!is_type(new_names, "character", self$n_vars())) {
                err_msg("set_names", "new_names", "the same length as # variants")
            }
            var_inds <- 0:(length(new_names) - 1)
            set_var_set_var_names(self$genomes, var_inds, new_names)
            invisible(self)
        },

        # Add one or more blank, named variants
        add_vars = function(new_names) {

            private$check_ptr()

            if (!is_type(new_names, "character")) {
                err_msg("add_vars", "new_names", "a character vector")
            }
            if (anyDuplicated(new_names) != 0) {
                err_msg("add_vars", "new_names", "a vector with no duplicates")
            }
            if (any(new_names %in% self$var_names())) {
                err_msg("add_vars", "new_names", "a vector containing no names",
                        "already present as variants names in the object being added to")
            }

            add_var_set_vars(self$genomes, new_names)

            invisible(self)
        },

        # Duplicate one or more variants by name
        dup_vars = function(var_names, new_names = NULL) {
            private$check_ptr()
            self_names <- self$var_names()
            if (!is_type(var_names, "character")) {
                err_msg("dup_vars", "var_names", "a character vector")
            }
            if (!all(var_names %in% self_names)) {
                err_msg("dup_vars", "var_names", "a character vector containing only",
                        "variant names present in the `variants` object.")
            }
            var_inds <- match(var_names, self_names) - 1

            if (!is.null(new_names) && !is_type(new_names, "character")) {
                err_msg("dup_vars", "new_names", "NULL or a character vector")
            }
            if (is.null(new_names)) {
                dup_n <- numeric(length(unique(var_names))) + 1
                names(dup_n) <- unique(var_names)
                new_names <- var_names
                for (i in 1:length(new_names)) {
                    new_names[i] <- paste0(var_names[i], "_dup", dup_n[[var_names[i]]])
                    dup_n[var_names[i]] <- dup_n[[var_names[i]]] + 1
                }
            }

            if (length(new_names) != length(var_names)) {
                err_msg("dup_vars", "new_names", "the same length as `var_names`")
            }
            if (any(new_names %in% self_names)) {
                err_msg("dup_seqs", "new_names", "a vector containing no elements",
                        "already present as a variant name in the variants object")
            }

            dup_var_set_vars(self$genomes, var_inds, new_names)

            invisible(self)
        },

        # Remove one or more variants by name
        rm_vars = function(var_names) {
            private$check_ptr()
            self_names <- self$var_names()
            if (!is_type(var_names, "character")) {
                err_msg("rm_vars", "var_names", "a character vector")
            }
            if (!all(var_names %in% self_names)) {
                err_msg("rm_vars", "var_names", "a vector of only names that are",
                        "present in the variants object")
            }
            if (anyDuplicated(var_names) != 0) {
                err_msg("rm_vars", "var_names", "a vector of non-duplicate names")
            }
            var_inds <- match(var_names, self_names) - 1
            remove_var_set_vars(self$genomes, var_inds)
            invisible(self)
        },


        # Mutations:
        add_sub = function(var_ind, seq_ind, pos, nt) {
            private$check_pos(var_ind, seq_ind, pos, "add_sub", "pos")
            if (!is_type(nt, "character", 1) || nchar(nt) != 1) {
                err_msg("add_sub", "nt", "a single character")
            }
            if (! nt %in% c("T", "C", "A", "G", "N")) {
                err_msg("add_sub", "nt", "one of \"T\", \"C\", \"A\", \"G\", or \"N\"")
            }
            add_substitution(self$genomes, var_ind - 1, seq_ind - 1, nt, pos - 1)
            invisible(self)
        },

        add_ins = function(var_ind, seq_ind, pos, nts) {
            private$check_pos(var_ind, seq_ind, pos, "add_ins", "pos")
            if (!is_type(nts, "character", 1)) {
                err_msg("add_ins", "nts", "a single string")
            }
            if (! all(strsplit(nts, "")[[1]] %in% c("T", "C", "A", "G", "N"))) {
                err_msg("add_ins", "nts", "string containing only \"T\", \"C\", \"A\",",
                        "\"G\", or \"N\"")
            }
            add_insertion(self$genomes, var_ind - 1, seq_ind - 1, nts, pos - 1)
            invisible(self)
        },

        add_del = function(var_ind, seq_ind, pos, n_nts) {
            private$check_pos(var_ind, seq_ind, pos, "add_del", "pos")
            if (!single_integer(n_nts, 1)) {
                err_msg("add_del", "n_nts", "a single integer >= 1")
            }
            add_deletion(self$genomes, var_ind - 1, seq_ind - 1, n_nts, pos - 1)
            invisible(self)
        }

    ),


    # __private__ ------
    private = list(
        # This should store a `XPtr<RefGenome>` to make sure it doesn't
        # go out of scope:
        reference = NULL,

        check_ptr = function() {
            if (!inherits(self$genomes, "externalptr")) {
                stop("\nSomehow the `genomes` field of this `variants` object ",
                     "has been converted to something other than an `externalptr` ",
                     "object. Re-make this `variants` object, don't edit the ",
                     "`genomes` field directly, and try again.", call. = FALSE)
            }
        },

        check_seq_ind = function(seq_ind, .fun_name) {
            private$check_ptr()
            if (!single_integer(seq_ind, 1, self$n_seqs())) {
                err_msg(.fun_name, "seq_ind", "integer in range [1, <# sequences>]")
            }
        },
        check_var_ind = function(var_ind, .fun_name) {
            private$check_ptr()
            if (!single_integer(var_ind, 1, self$n_vars())) {
                err_msg(.fun_name, "var_ind", "integer in range [1, <# variants>]")
            }
        },
        check_pos = function(var_ind, seq_ind, pos, .fun_name, .pos_name) {
            private$check_ptr()
            private$check_seq_ind(seq_ind, .fun_name)
            private$check_var_ind(var_ind, .fun_name)
            if (!single_integer(pos, 1, self$sizes(var_ind)[seq_ind])) {
                err_msg(.fun_name, .pos_name, "integer in range [1, <sequence size>]")
            }
        }
    )

)

variants$lock()





# >> CLASS mevo ----
#' An R6 class containing information needed for molecular evolution.
#'
#' This class is only used in `create_variants` to organize information.
#' It is not exported.
#'
#' @field Q A matrix of substitution rates for each nucleotide.
#' @field pi_tcag Vector of nucleotide equilibrium frequencies for "T", "C", "A", and
#'     "G", respectively.
#' @field insertion_rates Vector of insertion rates by length.
#' @field deletion_rates Vector of deletion rates by length.
#' @field gamma_mats List of matrices specifying "gamma distances" (see definition in
#'     `?create_mevo`) for each sequence.
#'
#' @section Methods:
#' \describe{
#'     \item{`mu()`}{Calculates the average overall mutation rate at equilibrium.}
#'     \item{`q()`}{Calculates the mutation rate for each nucleotide.}
#' }
#'
#' @return An object of class \code{mevo}.
#'
#' @docType class
#'
#' @seealso \code{\link{create_mevo}}
#'
#' @noRd
#'
#' @format An \code{\link[R6]{R6Class}} generator object
#'
#' @importFrom R6 R6Class
#'
mevo <- R6Class(

    "mevo",

    public = list(

        Q = NULL,
        pi_tcag = NULL,
        insertion_rates = NULL,
        deletion_rates = NULL,
        gamma_mats = NULL,
        region_size = NULL,

        initialize = function(sub_info,
                              insertion_rates,
                              deletion_rates,
                              gamma_mats,
                              region_size) {

            self$Q <- sub_info$Q
            self$pi_tcag <- sub_info$pi_tcag
            self$insertion_rates <- insertion_rates
            self$deletion_rates <- deletion_rates
            self$gamma_mats <- gamma_mats
            self$region_size <- region_size

        },


        print = function(digits = max(3, getOption("digits") - 3), ...) {
            fmt <- paste0("%.", digits, "f")
            cat("< Molecular evolution info >\n")

            cat("# Equilibrium densities:\n")
            cat("  ", sprintf(fmt, self$pi_tcag), "\n")

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
            prmatrix(self$Q, digits = digits,
                     rowlab = paste("  ", c("T", "C", "A", "G")),
                     collab = c("T", "C", "A", "G"))

            invisible(self)
        },


        # Average mutation rate
        mu = function() {
            # Indel rates (same for each nucleotide):
            indel <- sum(self$insertion_rates * 0.25) + sum(self$deletion_rates * 0.25)
            # Average mutation rate among all nucleotides:
            mu <- sum({rowSums(self$Q) + indel} * self$pi_tcag)
            return(mu)
        },

        # Overall mutation rate by nucleotide
        q = function() {
            # Indel rates (same for each nucleotide):
            indel <- sum(self$insertion_rates * 0.25) + sum(self$deletion_rates * 0.25)
            # Mutation rates by nucleotides:
            q <- rowSums(self$Q) + indel
            return(q)
        }

    ),


    private = list()

)

mevo$lock()




