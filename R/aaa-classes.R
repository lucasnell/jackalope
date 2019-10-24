# >> CLASS ref_genome----
#' An R6 class representing a reference genome.
#'
#'
#' This class should NEVER be created using `ref_genome$new`.
#' Only use `read_fasta` or `create_genome`.
#' This class wraps a pointer to a C++ object, which is why
#' there are no fields to manipulate directly.
#' All manipulations are done through this class's methods.
#'
#'
#'
#' @section Methods:
#' \strong{Viewing information:}
#' \describe{
#'     \item{`ptr()`}{View the pointer to the reference chromosome information.
#'         This is used internally by `jackalope` and shouldn't be of much use to users.}
#'     \item{`n_chroms()`}{View the number of chromosomes.}
#'     \item{`sizes()`}{View vector of chromosome sizes.}
#'     \item{`chrom_names()`}{View vector of chromosome names.}
#'     \item{`chrom(chrom_ind)`}{View a chromosome sequence string based on an index,
#'         `chrom_ind`.}
#'     \item{`gc_prop(chrom_ind, start, end)`}{View the GC proportion for a range within a
#'         reference chromosome.}
#'     \item{`nt_prop(nt, chrom_ind, start, end)`}{View the proportion of a range within a
#'         reference chromosome that is of nucleotide `nt`.}
#' }
#' \strong{Editing information:}
#' \describe{
#'     \item{`set_names(new_names)`}{Set names for all chromosomes.
#'         `new_names` is a character vector of what to change names to, and it must
#'         be the same length as the # chromosomes.}
#'     \item{`clean_names()`}{Clean chromosome names, converting `" :;=%,\\|/\"\'"`
#'         to `"_"`.}
#'     \item{`add_chroms(new_chroms, new_names = NULL)`}{Add one or more chromosomes
#'         directly. They can optionally be named (using `new_names`).
#'         Otherwise, their names are auto-generated.}
#'     \item{`rm_chroms(chrom_names)`}{Remove one or more chromosomes based on names in
#'         the `chrom_names` vector.}
#'     \item{`merge_chroms()`}{Merge all chromosomes into one after first shuffling
#'         their order.}
#'     \item{`filter_chroms(threshold, method)`}{Filter chromosomes by size
#'         (`method = "size"`) or for a proportion of total bases (`method = "prop"`).
#'         For the latter, chromosomes are first size-sorted, then the largest `N`
#'         chromosomes are retained that allow at least
#'         `threshold * sum(<all chromosome sizes>)` base pairs remaining after
#'         filtering.}
#'     \item{`replace_Ns(pi_tcag, n_threads = 1, show_progress = FALSE)`}{Replace
#'         `N`s in reference chromosome with nucleotides sampled with probabilities
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

        initialize = function(genome_ptr) {
            if (!inherits(genome_ptr, "externalptr")) {
                stop("\nWhen initializing a ref_genome object, you need to use ",
                     "an externalptr object.", call. = FALSE)
            }
            private$genome <- genome_ptr
        },

        print = function(...) {
            private$check_ptr()
            print_ref_genome(private$genome)
            invisible(self)
        },

        # ----------*
        # __view__ ----
        # ----------*
        # Pointer to underlying C++ object
        ptr = function() {
            private$check_ptr()
            return(private$genome)
        },

        # Get # chromosomes
        n_chroms = function() {
            private$check_ptr()
            return(view_ref_genome_nchroms(private$genome))
        },

        # Get vector of chromosome sizes
        sizes = function() {
            private$check_ptr()
            return(view_ref_genome_chrom_sizes(private$genome))
        },

        # Get vector of chromosome names
        chrom_names = function() {
            private$check_ptr()
            return(view_ref_genome_chrom_names(private$genome))
        },

        # Extract one reference chromosome
        chrom = function(chrom_ind) {
            private$check_ptr()
            if (!single_integer(chrom_ind, 1, self$n_chroms())) {
                err_msg("chrom", "chrom_ind", "integer in range [1, <# chromosomes>]")
            }
            return(view_ref_genome_chrom(private$genome, chrom_ind - 1))
        },
        # GC proportion for part of one reference chromosome
        gc_prop = function(chrom_ind, start, end) {
            private$check_pos(chrom_ind, start, "nt_prop", "start")
            private$check_pos(chrom_ind, end, "nt_prop", "end")
            if (end < start) err_msg("gc_prop", "end", ">= `start` arg")
            gcp <- view_ref_genome_gc_content(private$genome, chrom_ind - 1,
                                              start - 1, end - 1)
            return(gcp)
        },
        # Nucleotide content for part of one reference chromosome
        nt_prop = function(nt, chrom_ind, start, end) {
            private$check_pos(chrom_ind, start, "nt_prop", "start")
            private$check_pos(chrom_ind, end, "nt_prop", "end")
            if (end < start)err_msg("nt_prop", "end", ">= `start` arg")
            if (!is_type(nt, "character", 1) || nchar(nt) != 1) {
                err_msg("nt_prop", "nt", "a single character")
            }
            ntp <- view_ref_genome_nt_content(private$genome, nt, chrom_ind - 1,
                                              start - 1, end - 1)
            return(ntp)
        },

        # ----------*
        # __edit__ ----
        # ----------*
        # Change chromosome names
        set_names = function(new_names) {
            private$check_ptr()
            if (!is_type(new_names, "character", self$n_chroms())) {
                err_msg("set_names", "new_names", "the same length as # chromosomes")
            }
            chrom_inds <- 0:(length(new_names) - 1)
            set_ref_genome_chrom_names(private$genome, chrom_inds, new_names)
            invisible(self)
        },

        # Clean chromosome names, converting " :;=%,\\|/\"\'" to "_"
        clean_names = function() {
            private$check_ptr()
            clean_ref_genome_chrom_names(private$genome)
            invisible(self)
        },

        # Add one or more chromosomes
        add_chroms = function(new_chroms, new_names = NULL) {
            private$check_ptr()
            self_names <- self$chrom_names()
            if (!is_type(new_chroms, "character")) {
                err_msg("add_chroms", "new_chroms", "a character vector")
            }
            if (!is.null(new_names) && !is_type(new_names, "character")) {
                err_msg("add_chroms", "new_names", "NULL or a character vector")
            }
            if (is.null(new_names)) {
                N <- suppressWarnings(as.integer(gsub("chrom", "", self$chrom_names())))
                N <- if (any(!is.na(N))) 1 + max(N, na.rm = TRUE) else self$n_chroms()
                N <- N:(N+length(new_chroms)-1)
                new_names <- sprintf("chrom%i", N)
            }

            if (length(new_names) != length(new_chroms)) {
                err_msg("add_chroms","new_names", "the same length as `new_chroms`")
            }
            if (anyDuplicated(new_names) != 0) {
                err_msg("add_chroms","new_names", "a vector with no duplicates")
            }
            if (any(new_names %in% self$chrom_names())) {
                err_msg("add_chroms", "new_names", "a vector containing no elements",
                        "already present as a chromosome name in the reference genome")
            }

            add_ref_genome_chroms(private$genome, new_chroms, new_names)
            invisible(self)
        },

        # Remove one or more chromosomes by name
        rm_chroms = function(chrom_names) {
            private$check_ptr()
            self_names <- self$chrom_names()
            if (!is_type(chrom_names, "character")) {
                err_msg("rm_chroms", "chrom_names", "a character vector")
            }
            if (!all(chrom_names %in% self_names)) {
                err_msg("rm_chroms", "chrom_names", "a character vector of chromosome names",
                        "present in the `ref_genome` object. One or more of the names",
                        "provided weren't found")
            }
            if (anyDuplicated(chrom_names) != 0) {
                err_msg("rm_chroms", "chrom_names", "a vector of *non-duplicated* names")
            }
            chrom_inds <- match(chrom_names, self_names) - 1
            remove_ref_genome_chroms(private$genome, chrom_inds)
            invisible(self)
        },

        # Merge all ref_genome genome chromosomes into one
        merge_chroms = function() {
            private$check_ptr()
            merge_chromosomes_cpp(private$genome)
            invisible(self)
        },

        # Filter ref_genome chromosomes by size or for a proportion of total bases
        filter_chroms = function(threshold, method) {
            private$check_ptr()
            method <- match.arg(method, c("size", "prop"))
            min_chrom_size <- 0
            out_chrom_prop <- 0
            # Filling in the necessary parameter and checking for sensible inputs
            if (!single_number(threshold)) {
                err_msg("filter_chroms", "threshold", "a single number")
            }
            if (method == "size") {
                if (!single_integer(threshold, 1)) {
                    err_msg("filter_chroms", "threshold", "a single integer >= 1 if",
                            "filtering based on chromosome sizes")
                }
                min_chrom_size <- threshold
            } else {
                if (!single_number(threshold) || threshold >= 1 || threshold <= 0) {
                    err_msg("filter_chroms", "threshold", "a single number > 0 and < 1 if",
                            "filtering based on a proportion of total bases")
                }
                out_chrom_prop <- threshold
            }
            filter_chromosomes_cpp(private$genome, min_chrom_size, out_chrom_prop)
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

            replace_Ns_cpp(private$genome, pi_tcag, n_threads, show_progress)

            invisible(self)

        }


    ),

    # __private__ -----
    private = list(

        genome = NULL,

        check_ptr = function() {
            if (!inherits(private$genome, "externalptr")) {
                stop("\nSomehow the pointer to the underlying C++ object for this",
                     "`ref_genome` object has been converted to something other than ",
                     "an `externalptr` object---or was never set as one. ",
                     "You should re-make this `ref_genome` object using `read_fasta` ",
                     "or `create_genome`.", call. = FALSE)
            }
        },

        check_pos = function(chrom_ind, pos, .fun, .pos_name) {
            private$check_ptr()
            if (!single_integer(chrom_ind, 1, self$n_chroms())) {
                err_msg(.fun, "chrom_ind", "integer in range [1, <# chromosomes>]")
            }
            if (!single_integer(pos, 1, self$sizes()[chrom_ind])) {
                err_msg(.fun, .pos_name, "integer in range [1, <chromosome size>]")
            }
        }

    ),

    lock_class = TRUE

)









# >> CLASS variants----
#' An R6 class representing haploid variants from a reference genome.
#'
#' This class should NEVER be created using `variants$new`.
#' Only use `create_variants`.
#' This class wraps a pointer to a C++ object, which is why
#' there are no fields to manipulate directly.
#' All manipulations are done through this class's methods.
#'
#' @section Connections to `ref_genome` objects:
#' Regarding the `ref_genome` object you use to create a `variants` object, you should
#' note the following:
#'
#' \itemize{
#'     \item \strong{This point is the most important.}
#'         Both the `ref_genome` and `variants` objects use the same underlying
#'         C++ object to store reference genome information.
#'         Thus, if you make any changes to the `ref_genome` object, those changes will
#'         also show up in the `variants` object.
#'         For example, if you make a `variants` object named `V`
#'         based on an existing `ref_genome` object named `R`,
#'         then you merge chromosomes in `R`,
#'         `V` will now have merged chromosomes.
#'         If you've already started adding mutations to `V`,
#'         then all the indexes used to store those mutations will be inaccurate.
#'         So when you do anything with `V` later, your R session will crash
#'         or have errors.
#'         \strong{The lesson here is that you shouldn't edit the reference
#'         genome after using it to create variants.}
#'     \item If a `ref_genome` object is used to create a `variants`
#'         object, deleting the `ref_genome` object won't cause issues with
#'         the `variants` object.
#'         However, the `variants` class doesn't provide methods to edit
#'         chromosomes, so only remove the `ref_genome` object when you're done
#'         editing the reference genome.
#' }
#'
#' @section Methods:
#' \strong{Viewing information:}
#' \describe{
#'     \item{`ptr()`}{View the pointer to the variant information.
#'         This is used internally by `jackalope` and shouldn't be of much use to users.}
#'     \item{`n_chroms()`}{View the number of chromosomes.}
#'     \item{`n_vars()`}{View the number of variants.}
#'     \item{`sizes(var_ind)`}{View vector of chromosome sizes for a given variant.}
#'     \item{`chrom_names()`}{View vector of chromosome names.}
#'     \item{`var_names()`}{View vector of variant names.}
#'     \item{`chrom(var_ind, chrom_ind)`}{View a chromosome sequence string based on
#'         indices for the chromosome (`chrom_ind`) and variant (`var_ind`).}
#'     \item{`gc_prop(var_ind, chrom_ind, start, end)`}{View the GC proportion for a range
#'         within a variant chromosome.}
#'     \item{`nt_prop(nt, var_ind, chrom_ind, start, end)`}{View the proportion of a range
#'         within a variant chromosome that is of nucleotide `nt`.}
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
#'     \item{`add_sub(var_ind, chrom_ind, pos, nt)`}{Manually add a substitution
#'         for a given variant (`var_ind`), chromosome (`chrom_ind`), and position (`pos`).
#'         The reference nucleotide will be changed to `nt`, which should be a single
#'         character.}
#'     \item{`add_ins(var_ind, chrom_ind, pos, nts)`}{Manually add an insertion
#'         for a given variant (`var_ind`), chromosome (`chrom_ind`), and position (`pos`).
#'         The nucleotide(s) `nts` will be inserted after the designated position.}
#'     \item{`add_del(var_ind, chrom_ind, pos, n_nts)`}{Manually add a deletion
#'         for a given variant (`var_ind`), chromosome (`chrom_ind`), and position (`pos`).
#'         The designated number of nucleotides to delete (`n_nts`) will be deleted
#'         starting at `pos`, unless `pos` is near the chromosome end and doesn't have
#'         `n_nts` nucleotides to remove; it simply stops at the chromosome end in
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

        initialize = function(genomes_ptr, reference_ptr) {
            if (!inherits(genomes_ptr, "externalptr")) {
                stop("\nWhen initializing a variants object, you need to use ",
                     "an externalptr object.", call. = FALSE)
            }
            if (!inherits(reference_ptr, "externalptr")) {
                stop("\nWhen initializing a variants object, you need to use ",
                     "an externalptr object.", call. = FALSE)
            }
            private$genomes <- genomes_ptr
            private$reference <- reference_ptr
        },

        print = function() {
            private$check_ptr()
            print_var_set(private$genomes)
            invisible(self)
        },


        # ----------*
        # __view__ ----
        # ----------*
        # Pointer to underlying C++ object
        ptr = function() {
            private$check_ptr()
            return(private$genomes)
        },

        # Get # chromosomes
        n_chroms = function() {
            private$check_ptr()
            return(view_var_set_nchroms(private$genomes))
        },

        # Get # variants
        n_vars = function() {
            private$check_ptr()
            return(view_var_set_nvars(private$genomes))
        },

        # Get vector of chromosome sizes for one variant
        sizes = function(var_ind) {
            private$check_ptr()
            private$check_var_ind(var_ind, "sizes")
            return(view_var_genome_chrom_sizes(private$genomes, var_ind - 1))
        },

        # Get vector of chromosome names
        chrom_names = function() {
            stopifnot(inherits(private$reference, "externalptr"))
            return(view_ref_genome_chrom_names(private$reference))
        },

        # Get vector of variant names
        var_names = function() {
            private$check_ptr()
            return(view_var_set_var_names(private$genomes))
        },

        # Extract one variant chromosome
        chrom = function(var_ind, chrom_ind) {
            private$check_ptr()
            private$check_var_ind(var_ind, "chrom")
            private$check_chrom_ind(chrom_ind, "chrom")
            return(view_var_genome_chrom(private$genomes, var_ind - 1, chrom_ind - 1))
        },

        # GC proportion for part of one variant chromosome
        gc_prop = function(var_ind, chrom_ind, start, end) {
            private$check_pos(var_ind, chrom_ind, start, "gc_prop", "start")
            private$check_pos(var_ind, chrom_ind, end, "gc_prop", "end")
            if (end < start) err_msg("gc_prop", "end", ">= `start` arg")
            gcp <- view_var_set_gc_content(private$genomes, chrom_ind - 1, var_ind - 1,
                                           start - 1, end - 1)
            return(gcp)
        },

        # Nucleotide content for part of one reference chromosome
        nt_prop = function(nt, var_ind, chrom_ind, start, end) {
            private$check_pos(var_ind, chrom_ind, start, "nt_prop", "start")
            private$check_pos(var_ind, chrom_ind, end, "nt_prop", "end")
            if (end < start) err_msg("nt_prop", "end", ">= `start` arg")
            if (!is_type(nt, "character", 1) || nchar(nt) != 1) {
                err_msg("nt_prop", "nt", "a single character")
            }
            ntp <- view_var_set_nt_content(private$genomes, nt, chrom_ind - 1, var_ind - 1,
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
            set_var_set_var_names(private$genomes, var_inds, new_names)
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

            add_var_set_vars(private$genomes, new_names)

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
                err_msg("dup_vars", "new_names", "a vector containing no elements",
                        "already present as a variant name in the variants object")
            }

            dup_var_set_vars(private$genomes, var_inds, new_names)

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
            remove_var_set_vars(private$genomes, var_inds)
            invisible(self)
        },


        # Mutations:
        add_sub = function(var_ind, chrom_ind, pos, nt) {
            private$check_pos(var_ind, chrom_ind, pos, "add_sub", "pos")
            if (!is_type(nt, "character", 1) || nchar(nt) != 1) {
                err_msg("add_sub", "nt", "a single character")
            }
            if (! nt %in% c("T", "C", "A", "G", "N")) {
                err_msg("add_sub", "nt", "one of \"T\", \"C\", \"A\", \"G\", or \"N\"")
            }
            add_substitution(private$genomes, var_ind - 1, chrom_ind - 1, nt, pos - 1)
            invisible(self)
        },

        add_ins = function(var_ind, chrom_ind, pos, nts) {
            private$check_pos(var_ind, chrom_ind, pos, "add_ins", "pos")
            if (!is_type(nts, "character", 1)) {
                err_msg("add_ins", "nts", "a single string")
            }
            if (! all(strsplit(nts, "")[[1]] %in% c("T", "C", "A", "G", "N"))) {
                err_msg("add_ins", "nts", "string containing only \"T\", \"C\", \"A\",",
                        "\"G\", or \"N\"")
            }
            add_insertion(private$genomes, var_ind - 1, chrom_ind - 1, nts, pos - 1)
            invisible(self)
        },

        add_del = function(var_ind, chrom_ind, pos, n_nts) {
            private$check_pos(var_ind, chrom_ind, pos, "add_del", "pos")
            if (!single_integer(n_nts, 1)) {
                err_msg("add_del", "n_nts", "a single integer >= 1")
            }
            add_deletion(private$genomes, var_ind - 1, chrom_ind - 1, n_nts, pos - 1)
            invisible(self)
        }

    ),


    # __private__ ------
    private = list(

        genomes = NULL,

        # This should store a `XPtr<RefGenome>` to make sure it doesn't
        # go out of scope:
        reference = NULL,

        check_ptr = function() {
            if (!inherits(private$genomes, "externalptr")) {
                stop("\nSomehow the pointer to the underlying C++ object for this ",
                     "`variants` object has been converted to something other than ",
                     "an `externalptr` object---or was never set as one. ",
                     "You should re-make this `variants` object using ",
                     "`create_variants`.", call. = FALSE)
            }
        },

        check_chrom_ind = function(chrom_ind, .fun_name) {
            private$check_ptr()
            if (!single_integer(chrom_ind, 1, self$n_chroms())) {
                err_msg(.fun_name, "chrom_ind", "integer in range [1, <# chromosomes>]")
            }
        },
        check_var_ind = function(var_ind, .fun_name) {
            private$check_ptr()
            if (!single_integer(var_ind, 1, self$n_vars())) {
                err_msg(.fun_name, "var_ind", "integer in range [1, <# variants>]")
            }
        },
        check_pos = function(var_ind, chrom_ind, pos, .fun_name, .pos_name) {
            private$check_ptr()
            private$check_chrom_ind(chrom_ind, .fun_name)
            private$check_var_ind(var_ind, .fun_name)
            if (!single_integer(pos, 1, self$sizes(var_ind)[chrom_ind])) {
                err_msg(.fun_name, .pos_name, "integer in range [1, <chromosome size>]")
            }
        }
    ),

    lock_class = TRUE

)




