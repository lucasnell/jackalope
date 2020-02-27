# >> CLASS ref_genome----
#' R6 Class Representing a Reference Genome
#'
#' @description
#' Interactive wrapper for a pointer to a C++ object that stores reference genome
#' information.
#'
#'
#' @details
#' This class should NEVER be created using `ref_genome$new`.
#' Only use `read_fasta` or `create_genome`.
#' Because this class wraps a pointer to a C++ object, there are no fields to
#' manipulate directly.
#' All manipulations are done through this class's methods.
#'
#'
#' @seealso \code{\link{read_fasta}} \code{\link{create_genome}}
#'
#' @export
#'
#' @param chrom_ind Index for the focal chromosome.
#' @param start Point on the chromosome at which to start the calculation
#'     (inclusive).
#' @param end Point on the chromosome at which to end the calculation
#'     (inclusive).
#'
#' @importFrom R6 R6Class
#'
ref_genome <- R6Class(

    "ref_genome",

    public = list(


        #' @description
        #' Do NOT use this; only use `read_fasta` or `create_genome` to make a
        #' new `ref_genome`.
        #'
        #' @param genome_ptr An `externalptr` object pointing to a C++ object that stores
        #'     the information about the reference genome.
        #'
        initialize = function(genome_ptr) {
            if (!inherits(genome_ptr, "externalptr")) {
                stop("\nWhen initializing a ref_genome object, you need to use ",
                     "an externalptr object.", call. = FALSE)
            }
            private$genome <- genome_ptr
        },

        #' @description
        #' Print a `ref_genome` object.
        print = function() {
            private$check_ptr()
            print_ref_genome(private$genome)
            invisible(self)
        },

        # ----------*
        # __view__ ----
        # ----------*

        #' @description
        #' View pointer to underlying C++ object (this is not useful to end users).
        #'
        #' @return An `externalptr` object.
        #'
        ptr = function() {
            private$check_ptr()
            return(private$genome)
        },

        #' @description
        #' View number of chromosomes.
        #'
        #' @return Integer number of chromosomes.
        #'
        n_chroms = function() {
            private$check_ptr()
            return(view_ref_genome_nchroms(private$genome))
        },

        #' @description
        #' View chromosome sizes.
        #'
        #' @return Integer vector of chromosome sizes.
        #'
        sizes = function() {
            private$check_ptr()
            return(view_ref_genome_chrom_sizes(private$genome))
        },

        #' @description
        #' View chromosome names.
        #'
        #' @return Character vector of chromosome names.
        #'
        chrom_names = function() {
            private$check_ptr()
            return(view_ref_genome_chrom_names(private$genome))
        },

        #' @description
        #' View one reference chromosome.
        #'
        #' @return A single string representing the chosen chromosome's DNA sequence.
        #'
        chrom = function(chrom_ind) {
            private$check_ptr()
            if (!single_integer(chrom_ind, 1, self$n_chroms())) {
                err_msg("chrom", "chrom_ind", "integer in range [1, <# chromosomes>]")
            }
            return(view_ref_genome_chrom(private$genome, chrom_ind - 1))
        },


        #' @description
        #' View GC proportion for part of one reference chromosome.
        #'
        #'
        #' @return A double in the range `[0,1]` representing the proportion of DNA
        #' sequence that is either `G` or `C`.
        #'
        gc_prop = function(chrom_ind, start, end) {
            private$check_pos(chrom_ind, start, "nt_prop", "start")
            private$check_pos(chrom_ind, end, "nt_prop", "end")
            if (end < start) err_msg("gc_prop", "end", ">= `start` arg")
            gcp <- view_ref_genome_gc_content(private$genome, chrom_ind - 1,
                                              start - 1, end - 1)
            return(gcp)
        },

        #' @description
        #' View nucleotide content for part of one reference chromosome
        #'
        #' @param nt Which nucleotide to calculate the proportion that the DNA
        #'     sequence is made of. Must be one of `T`, `C`, `A`, `G`, or `N`.
        #'
        #' @return A double in the range `[0,1]` representing the proportion of DNA
        #' sequence that is `nt`.
        #'
        nt_prop = function(nt, chrom_ind, start, end) {
            private$check_pos(chrom_ind, start, "nt_prop", "start")
            private$check_pos(chrom_ind, end, "nt_prop", "end")
            if (end < start) err_msg("nt_prop", "end", ">= `start` arg")
            if (!is_type(nt, "character", 1) || nchar(nt) != 1) {
                err_msg("nt_prop", "nt", "a single character")
            }
            if (!nt %in% c("T", "C", "A", "G", "N")) {
                err_msg("nt_prop", "nt", 'either "T", "C", "A", "G", or "N"')
            }
            ntp <- view_ref_genome_nt_content(private$genome, nt, chrom_ind - 1,
                                              start - 1, end - 1)
            return(ntp)
        },


        # ----------*
        # __edit__ ----
        # ----------*


        #' @description
        #' Change chromosome names.
        #'
        #' @param new_names Vector of new names to use. This must be the same length as
        #'     the number of current names.
        #'
        #' @return This `R6` object, invisibly.
        #'
        #' @examples
        #' ref <- create_genome(4, 10)
        #' ref$set_names(c("a", "b", "c", "d"))
        #'
        set_names = function(new_names) {
            private$check_ptr()
            if (!is_type(new_names, "character", self$n_chroms())) {
                err_msg("set_names", "new_names", "the same length as # chromosomes")
            }
            chrom_inds <- 0:(length(new_names) - 1)
            set_ref_genome_chrom_names(private$genome, chrom_inds, new_names)
            invisible(self)
        },

        #' @description
        #' Clean chromosome names, converting `" :;=%,\\|/\"\'"` to `"_"`.
        #'
        #' @return This `R6` object, invisibly.
        #'
        #' @examples
        #' ref <- create_genome(4, 10)
        #' ref$set_names(c("a:", "b|", "c;", "d'"))
        #' ref$clean_names()
        #'
        clean_names = function() {
            private$check_ptr()
            clean_ref_genome_chrom_names(private$genome)
            invisible(self)
        },

        #' @description
        #' Add one or more chromosomes.
        #'
        #' @param new_chroms Character vector of DNA strings representing new chromosomes.
        #' @param new_names Optional character vector of names for the new chromosomes.
        #'     It should be the same length as `new_chroms`.
        #'     If `NULL`, new names will be automatically generated. Defaults to `NULL`.
        #'
        #' @return This `R6` object, invisibly.
        #'
        #' @examples
        #' ref <- create_genome(4, 10)
        #' ref$add_chroms("TCAGTCAG")
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

        #' @description
        #' Remove one or more chromosomes by name
        #'
        #' @param chrom_names Vector of the name(s) of the chromosome(s) to remove.
        #'
        #' @return This `R6` object, invisibly.
        #'
        #' @examples
        #' ref <- create_genome(4, 10)
        #' ref$set_names(c("a", "b", "c", "d"))
        #' ref$rm_chroms("b")
        #'
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


        #' @description
        #' Merge chromosomes into one.
        #'
        #' @param chrom_names Vector of the names of the chromosomes to merge into one.
        #'     Duplicates are not allowed, and chromosomes are merged in the order
        #'     they're provided.
        #'     If this is `NULL`, then all chromosomes are merged after first
        #'     shuffling their order.
        #'
        #' @return This `R6` object, invisibly.
        #'
        #' @examples
        #' ref <- create_genome(4, 10)
        #' ref$merge_chroms(ref$chrom_names()[1:2])
        #' ref$merge_chroms(NULL)
        #'
        merge_chroms = function(chrom_names) {
            private$check_ptr()
            if (!is.null(chrom_names) && !is_type(chrom_names, "character")) {
                err_msg("merge_chroms", "chrom_names", "NULL or a character vector")
            }
            if (is.null(chrom_names)) {
                merge_all_chromosomes_cpp(private$genome)
            } else {
                if (length(chrom_names) < 2) {
                    err_msg("rm_chroms", "chrom_names", "a vector of >1 names")
                }
                self_names <- self$chrom_names()
                if (!all(chrom_names %in% self_names)) {
                    err_msg("merge_chroms", "chrom_names",
                            "a character vector of chromosome names",
                            "present in the `ref_genome` object. One or more of the names",
                            "provided weren't found")
                }
                if (anyDuplicated(chrom_names) != 0) {
                    err_msg("rm_chroms", "chrom_names",
                            "a vector of *non-duplicated* names")
                }
                chrom_inds <- match(chrom_names, self_names) - 1
                merge_chromosomes_cpp(private$genome, chrom_inds)
            }
            invisible(self)
        },


        #' @description
        #' Filter chromosomes by size or for a proportion of total bases.
        #'
        #' @param threshold Number used as a threshold. If `method == "size"`,
        #'     then this is the minimum length of a chromosome that will remain after
        #'     filtering.
        #'     If `method == "prop"`, chromosomes are first size-sorted, then
        #'     the largest `N` chromosomes are retained that allow at least
        #'     `threshold * sum(<all chromosome sizes>)` base pairs remaining after
        #'     filtering.
        #' @param method String indicating which filter method to use: chromosome size
        #'     (`method = "size"`) or proportion of total bases (`method = "prop"`).
        #'
        #' @return This `R6` object, invisibly.
        #'
        #' @examples
        #' ref <- create_genome(4, 100, 50)
        #' ref$filter_chroms(90, "size")
        #' ref$filter_chroms(0.4, "prop")
        #'
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

        #' @description
        #' Replace `N`s in the reference genome.
        #'
        #' @param pi_tcag Numeric vector (length 4) indicating the sampling weights
        #'     for `T`, `C`, `A`, and `G`, respectively, for generating new nucleotides
        #'     with which to replace the `N`s.
        #' @param n_threads Optional integer specifying the threads to use.
        #'     Ignored if the package wasn't compiled with OpenMP. Defaults to `1`.
        #' @param show_progress Optional logical indicating whether to show a
        #'     progress bar. Defaults to `FALSE`.
        #'
        #' @return This `R6` object, invisibly.
        #'
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

    lock_class = TRUE,
    cloneable = FALSE

)






# >> CLASS haplotypes----
#' An R6 Class Representing Haploid Variants
#'
#' @description
#' Interactive wrapper for a pointer to a C++ object that stores information about
#' variant haplotypes from a single reference genome.
#'
#' @details
#' This class should NEVER be created using `haplotypes$new`.
#' Only use `create_haplotypes`.
#' Because this class wraps a pointer to a C++ object, there are no fields to
#' manipulate directly.
#' All manipulations are done through this class's methods.
#'
#' @section Connections to `ref_genome` objects:
#' Regarding the `ref_genome` object you use to create a `haplotypes` object, you should
#' note the following:
#'
#' \itemize{
#'     \item \strong{This point is the most important.}
#'         Both the `ref_genome` and `haplotypes` objects use the same underlying
#'         C++ object to store reference genome information.
#'         Thus, if you make any changes to the `ref_genome` object, those changes will
#'         also show up in the `haplotypes` object.
#'         For example, if you make a `haplotypes` object named `V`
#'         based on an existing `ref_genome` object named `R`,
#'         then you merge chromosomes in `R`,
#'         `V` will now have merged chromosomes.
#'         If you've already started adding mutations to `V`,
#'         then all the indexes used to store those mutations will be inaccurate.
#'         So when you do anything with `V` later, your R session will crash
#'         or have errors.
#'         \strong{The lesson here is that you shouldn't edit the reference
#'         genome after using it to create haplotypes.}
#'     \item If a `ref_genome` object is used to create a `haplotypes`
#'         object, deleting the `ref_genome` object won't cause issues with
#'         the `haplotypes` object.
#'         However, the `haplotypes` class doesn't provide methods to edit
#'         chromosomes, so only remove the `ref_genome` object when you're done
#'         editing the reference genome.
#' }
#'
#'
#' @seealso \code{\link{create_haplotypes}}
#'
#' @export
#'
#'
#' @param chrom_ind Index for the focal chromosome.
#' @param hap_ind Index for the focal haplotype.
#' @param start Point on the chromosome at which to start the calculation
#'     (inclusive).
#' @param end Point on the chromosome at which to end the calculation
#'     (inclusive).
#' @param pos Position at which to add the mutation.
#'
#'
#' @importFrom R6 R6Class
#'
#'
haplotypes <- R6Class(

    "haplotypes",

    public = list(


        #' @description
        #' Do NOT use this; only use `create_haplotypes` to make new `haplotypes`.
        #'
        #' @param genomes_ptr An `externalptr` object pointing to a C++ object that
        #'     stores the information about the haplotypes.
        #' @param reference_ptr An `externalptr` object pointing to a C++ object that
        #'     stores the information about the reference genome.
        #'
        initialize = function(genomes_ptr, reference_ptr) {
            if (!inherits(genomes_ptr, "externalptr")) {
                stop("\nWhen initializing a haplotypes object, you need to use ",
                     "an externalptr object.", call. = FALSE)
            }
            if (!inherits(reference_ptr, "externalptr")) {
                stop("\nWhen initializing a haplotypes object, you need to use ",
                     "an externalptr object.", call. = FALSE)
            }
            private$genomes <- genomes_ptr
            private$reference <- reference_ptr
        },

        #' @description
        #' Print a `haplotypes` object.
        print = function() {
            private$check_ptr()
            print_hap_set(private$genomes)
            invisible(self)
        },


        # ----------*
        # __view__ ----
        # ----------*

        #' @description
        #' View pointer to underlying C++ object (this is not useful to end users).
        #'
        #' @return An `externalptr` object.
        #'
        ptr = function() {
            private$check_ptr()
            return(private$genomes)
        },

        #' @description
        #' View number of chromosomes.
        #'
        #' @return Integer number of chromosomes.
        #'
        n_chroms = function() {
            private$check_ptr()
            return(view_hap_set_nchroms(private$genomes))
        },

        #' @description
        #' View number of haplotypes.
        #'
        #' @return Integer number of haplotypes.
        #'
        n_haps = function() {
            private$check_ptr()
            return(view_hap_set_nhaps(private$genomes))
        },


        #' @description
        #' View chromosome sizes for one haplotype.
        #'
        #' @return Integer vector of chromosome sizes for focal haplotype.
        #'
        sizes = function(hap_ind) {
            private$check_ptr()
            private$check_hap_ind(hap_ind, "sizes")
            return(view_hap_genome_chrom_sizes(private$genomes, hap_ind - 1))
        },

        #' @description
        #' View chromosome names.
        #'
        #' @return Character vector of chromosome names.
        #'
        chrom_names = function() {
            stopifnot(inherits(private$reference, "externalptr"))
            return(view_ref_genome_chrom_names(private$reference))
        },

        #' @description
        #' View haplotype names.
        #'
        #' @return Character vector of haplotype names.
        #'
        hap_names = function() {
            private$check_ptr()
            return(view_hap_set_hap_names(private$genomes))
        },

        #' @description
        #' View one haplotype chromosome.
        #'
        #' @return A single string representing the chosen haplotype chromosome's DNA
        #' sequence.
        #'
        chrom = function(hap_ind, chrom_ind) {
            private$check_ptr()
            private$check_hap_ind(hap_ind, "chrom")
            private$check_chrom_ind(chrom_ind, "chrom")
            return(view_hap_genome_chrom(private$genomes, hap_ind - 1, chrom_ind - 1))
        },


        #' @description
        #' View GC proportion for part of one haplotype chromosome.
        #'
        #' @return A double in the range `[0,1]` representing the proportion of DNA
        #' sequence that is either `G` or `C`.
        #'
        gc_prop = function(hap_ind, chrom_ind, start, end) {
            private$check_pos(hap_ind, chrom_ind, start, "gc_prop", "start")
            private$check_pos(hap_ind, chrom_ind, end, "gc_prop", "end")
            if (end < start) err_msg("gc_prop", "end", ">= `start` arg")
            gcp <- view_hap_set_gc_content(private$genomes, chrom_ind - 1, hap_ind - 1,
                                           start - 1, end - 1)
            return(gcp)
        },

        #' @description
        #' View nucleotide content for part of one haplotype chromosome
        #'
        #' @param nt Which nucleotide to calculate the proportion that the DNA
        #'     sequence is made of. Must be one of `T`, `C`, `A`, `G`, or `N`.
        #'
        #' @return A double in the range `[0,1]` representing the proportion of DNA
        #' sequence that is `nt`.
        #'
        nt_prop = function(nt, hap_ind, chrom_ind, start, end) {
            private$check_pos(hap_ind, chrom_ind, start, "nt_prop", "start")
            private$check_pos(hap_ind, chrom_ind, end, "nt_prop", "end")
            if (end < start) err_msg("nt_prop", "end", ">= `start` arg")
            if (!is_type(nt, "character", 1) || nchar(nt) != 1) {
                err_msg("nt_prop", "nt", "a single character")
            }
            if (!nt %in% c("T", "C", "A", "G", "N")) {
                err_msg("nt_prop", "nt", 'either "T", "C", "A", "G", or "N"')
            }
            ntp <- view_hap_set_nt_content(private$genomes, nt, chrom_ind - 1, hap_ind - 1,
                                           start - 1, end - 1)
            return(ntp)
        },


        # ----------*
        # __edit__ ----
        # ----------*

        #' @description
        #' Change haplotype names.
        #'
        #' @param new_names Vector of new names to use. This must be the same length as
        #'     the number of current names.
        #'
        #' @return This `R6` object, invisibly.
        #'
        set_names = function(new_names) {
            private$check_ptr()
            if (!is_type(new_names, "character", self$n_haps())) {
                err_msg("set_names", "new_names", "the same length as # haplotypes")
            }
            hap_inds <- 0:(length(new_names) - 1)
            set_hap_set_hap_names(private$genomes, hap_inds, new_names)
            invisible(self)
        },

        #' @description
        #' Add one or more blank, named haplotypes
        #'
        #' @param new_names Vector of name(s) for the new haplotype(s).
        #'
        #' @return This `R6` object, invisibly.
        #'
        add_haps = function(new_names) {

            private$check_ptr()

            if (!is_type(new_names, "character")) {
                err_msg("add_haps", "new_names", "a character vector")
            }
            if (anyDuplicated(new_names) != 0) {
                err_msg("add_haps", "new_names", "a vector with no duplicates")
            }
            if (any(new_names %in% self$hap_names())) {
                err_msg("add_haps", "new_names", "a vector containing no names",
                        "already present as haplotypes names in the object being added to")
            }

            add_hap_set_haps(private$genomes, new_names)

            invisible(self)
        },


        #' @description
        #' Duplicate one or more haplotypes by name.
        #'
        #' @param hap_names Vector of existing haplotype name(s) that you want to
        #'     duplicate.
        #' @param new_names Optional vector specifying the names of the duplicates.
        #'     If `NULL`, their names are auto-generated. Defaults to `NULL`.
        #'
        #' @return This `R6` object, invisibly.
        #'
        dup_haps = function(hap_names, new_names = NULL) {
            private$check_ptr()
            self_names <- self$hap_names()
            if (!is_type(hap_names, "character")) {
                err_msg("dup_haps", "hap_names", "a character vector")
            }
            if (!all(hap_names %in% self_names)) {
                err_msg("dup_haps", "hap_names", "a character vector containing only",
                        "haplotype names present in the `haplotypes` object.")
            }
            hap_inds <- match(hap_names, self_names) - 1

            if (!is.null(new_names) && !is_type(new_names, "character")) {
                err_msg("dup_haps", "new_names", "NULL or a character vector")
            }
            if (is.null(new_names)) {
                dup_n <- numeric(length(unique(hap_names))) + 1
                names(dup_n) <- unique(hap_names)
                new_names <- hap_names
                for (i in 1:length(new_names)) {
                    new_names[i] <- paste0(hap_names[i], "_dup", dup_n[[hap_names[i]]])
                    dup_n[hap_names[i]] <- dup_n[[hap_names[i]]] + 1
                }
            }

            if (length(new_names) != length(hap_names)) {
                err_msg("dup_haps", "new_names", "the same length as `hap_names`")
            }
            if (any(new_names %in% self_names)) {
                err_msg("dup_haps", "new_names", "a vector containing no elements",
                        "already present as a haplotype name in the haplotypes object")
            }

            dup_hap_set_haps(private$genomes, hap_inds, new_names)

            invisible(self)
        },

        #' @description
        #' Remove one or more haplotypes by name.
        #'
        #' @param hap_names Vector of existing haplotype name(s) that you want to remove.
        #'
        #' @return This `R6` object, invisibly.
        #'
        rm_haps = function(hap_names) {
            private$check_ptr()
            self_names <- self$hap_names()
            if (!is_type(hap_names, "character")) {
                err_msg("rm_haps", "hap_names", "a character vector")
            }
            if (!all(hap_names %in% self_names)) {
                err_msg("rm_haps", "hap_names", "a vector of only names that are",
                        "present in the haplotypes object")
            }
            if (anyDuplicated(hap_names) != 0) {
                err_msg("rm_haps", "hap_names", "a vector of non-duplicate names")
            }
            hap_inds <- match(hap_names, self_names) - 1
            remove_hap_set_haps(private$genomes, hap_inds)
            invisible(self)
        },

        #' @description
        #' Manually add a substitution.
        #'
        #' @param nt Single character representing the nucleotide to change the
        #'     current one to.
        #'
        #' @return This `R6` object, invisibly.
        #'
        add_sub = function(hap_ind, chrom_ind, pos, nt) {
            private$check_pos(hap_ind, chrom_ind, pos, "add_sub", "pos")
            if (!is_type(nt, "character", 1) || nchar(nt) != 1) {
                err_msg("add_sub", "nt", "a single character")
            }
            if (! nt %in% c("T", "C", "A", "G", "N")) {
                err_msg("add_sub", "nt", "one of \"T\", \"C\", \"A\", \"G\", or \"N\"")
            }
            add_substitution(private$genomes, hap_ind - 1, chrom_ind - 1, nt, pos - 1)
            invisible(self)
        },

        #' @description
        #' Manually add an insertion.
        #'
        #' @param nts String representing the nucleotide(s) that will be inserted after
        #'     the designated position.
        #'
        #' @return This `R6` object, invisibly.
        #'
        add_ins = function(hap_ind, chrom_ind, pos, nts) {
            private$check_pos(hap_ind, chrom_ind, pos, "add_ins", "pos")
            if (!is_type(nts, "character", 1)) {
                err_msg("add_ins", "nts", "a single string")
            }
            if (! all(strsplit(nts, "")[[1]] %in% c("T", "C", "A", "G", "N"))) {
                err_msg("add_ins", "nts", "string containing only \"T\", \"C\", \"A\",",
                        "\"G\", or \"N\"")
            }
            add_insertion(private$genomes, hap_ind - 1, chrom_ind - 1, nts, pos - 1)
            invisible(self)
        },


#'     \item{`add_del(hap_ind, chrom_ind, pos, n_nts)`}{Manually add a deletion
#'         for a given haplotype (`hap_ind`), chromosome (`chrom_ind`), and position (`pos`).
#'         The designated number of nucleotides to delete (`n_nts`) will be deleted
#'         starting at `pos`, unless `pos` is near the chromosome end and doesn't have
#'         `n_nts` nucleotides to remove; it simply stops at the chromosome end in
#'         this case.}


        #' @description
        #' Manually add a deletion.
        #'
        #' @param n_nts Single integer specifying the number of nucleotides to delete.
        #'     These will be deleted starting at `pos`.
        #'     If `pos` is near the chromosome end and doesn't have `n_nts` nucleotides
        #'     to remove, it simply removes nucleotides from `pos` to the chromosome end.
        #'
        #' @return This `R6` object, invisibly.
        #'
        add_del = function(hap_ind, chrom_ind, pos, n_nts) {
            private$check_pos(hap_ind, chrom_ind, pos, "add_del", "pos")
            if (!single_integer(n_nts, 1)) {
                err_msg("add_del", "n_nts", "a single integer >= 1")
            }
            add_deletion(private$genomes, hap_ind - 1, chrom_ind - 1, n_nts, pos - 1)
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
                     "`haplotypes` object has been converted to something other than ",
                     "an `externalptr` object---or was never set as one. ",
                     "You should re-make this `haplotypes` object using ",
                     "`create_haplotypes`.", call. = FALSE)
            }
        },

        check_chrom_ind = function(chrom_ind, .fun_name) {
            private$check_ptr()
            if (!single_integer(chrom_ind, 1, self$n_chroms())) {
                err_msg(.fun_name, "chrom_ind", "integer in range [1, <# chromosomes>]")
            }
        },
        check_hap_ind = function(hap_ind, .fun_name) {
            private$check_ptr()
            if (!single_integer(hap_ind, 1, self$n_haps())) {
                err_msg(.fun_name, "hap_ind", "integer in range [1, <# haplotypes>]")
            }
        },
        check_pos = function(hap_ind, chrom_ind, pos, .fun_name, .pos_name) {
            private$check_ptr()
            private$check_chrom_ind(chrom_ind, .fun_name)
            private$check_hap_ind(hap_ind, .fun_name)
            if (!single_integer(pos, 1, self$sizes(hap_ind)[chrom_ind])) {
                err_msg(.fun_name, .pos_name, "integer in range [1, <chromosome size>]")
            }
        }
    ),

    lock_class = TRUE,
    cloneable = FALSE

)




