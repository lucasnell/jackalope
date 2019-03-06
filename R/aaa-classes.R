# >> ref_genome class----
#' An R6 class representing a reference genome.
#'
#'
#' \emph{Note:} Do NOT change fields in this class directly.
#' It will cause your R session to do bad things.
#' (Ever seen the bomb popup on RStudio? Manually mess with these fields and you
#' surely will.)
#' For safe ways of manipulating the reference genome, see the "Methods" section.
#'
#'
#' @field genome An \code{externalptr} to a C++ object storing the sequences
#'     representing the genome.
#' @field digests An \code{externalptr} to a C++ object storing the digestion of the
#'     genome, if a digestion has been carried out. It's \code{NULL} otherwise.
#'
#' @section Methods:
#' \describe{
#'     \item{`n_seqs()`}{View the number of sequences.}
#'     \item{`sizes()`}{View vector of sequence sizes.}
#'     \item{`names()`}{View vector of sequence names.}
#'     \item{`extract_seq(seq_ind)`}{Extract a sequence string based on an index,
#'         `seq_ind`.}
#'     \item{`set_names(names)`}{Set names for all sequences.}
#'     \item{`rm_seqs(seq_names)`}{Remove one or more sequences based on names in
#'         the `seq_names` vector.}
#'     \item{`merge_seqs()`}{Merge all sequences into one.}
#'     \item{`filter_seqs(threshold, method)`}{Filter sequences by size
#'         (`method = "size"`) or for a proportion of total bases `method = "prop"`.
#'         For the latter, sequences are first size-sorted, then the largest `N`
#'         sequences are retained that allow at least
#'         `threshold * sum(<all sequence sizes>)` base pairs remaining after
#'         filtering.}
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
        extract_seq = function(seq_ind) {
            private$check_ptr()
            if (seq_ind > self$n_seqs() | seq_ind < 1) {
                stop("seq_ind arg must be in range [1, <# sequences>]", call. = FALSE)
            }
            return(view_ref_genome_seq(self$genome, seq_ind - 1))
        },

        # ----------*
        # __edit__ ----
        # ----------*
        # Change sequence names
        set_names = function(names) {
            private$check_ptr()
            if (length(names) != self$n_seqs()) {
                stop("names arg must be the same length as # sequences", call. = FALSE)
            }
            seq_inds <- 0:(length(names) - 1)
            set_ref_genome_seq_names(self$genome, seq_inds, names)
            invisible(self)
        },

        # Remove one or more sequences by name
        rm_seqs = function(seq_names) {
            private$check_ptr()
            self_names <- self$names()
            if (!all(seq_names %in% self_names)) {
                stop("not all provided seq_names are in this genome.", call. = FALSE)
            }
            if (anyDuplicated(seq_names) != 0) {
                stop("one or more seq_names are duplicate.", call. = FALSE)
            }
            seq_inds <- match(seq_names, self_names) - 1
            remove_ref_genome_seqs(self$genome, seq_inds)
            invisible(self)
        },

        # Merge all ref_genome genome sequences into one
        merge_seqs = function() {
            private$check_ptr()
            merge_sequences(self$genome)
            invisible(self)
        },

        # Filter ref_genome sequences by size or for a proportion of total bases
        filter_seqs = function(threshold, method) {
            private$check_ptr()
            method <- match.arg(method, c("size", "prop"))
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

    private = list(

        check_ptr = function() stopifnot(inherits(self$genome, "externalptr"))

    )

)



ref_genome$lock()




# >> mevo class----
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
#' @section Methods:
#' \describe{
#'     \item{`mu()`}{Calculates the average overall mutation rate at equilibrium.}
#'     \item{`q()`}{Calculates the mutation rate for each nucleotide.}
#'     \item{`to_ptr()`}{Converts information in this object to a C++ pointer.
#'         You shouldn't need to use this. Ever.}
#' }
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
            dim(self$gamma_mats) <- NULL  # to make it a list instead of matrix
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
        }

    ),


    private = list()

)

mevo$lock()






# >> variants class----
#' An R6 class representing haploid variants from a reference genome.
#'
#' \emph{Note:} Do NOT change fields in this class directly. It will cause your
#' R session to do bad things.
#' (Ever seen the bomb popup on RStudio? Manually mess with these fields and you
#' surely will.)
#' For safe ways of manipulating the variants' information, see the "Methods" section.
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
#' @section Methods:
#' \describe{
#'     \item{`n_seqs()`}{View the number of sequences.}
#'     \item{`n_vars()`}{View the number of variants.}
#'     \item{`sizes(var_ind)`}{View vector of sequence sizes for a given variant.}
#'     \item{`seq_names()`}{View vector of sequence names.}
#'     \item{`var_names()`}{View vector of variant names.}
#'     \item{`extract_seq(var_ind, seq_ind)`}{Extract a sequence string based on
#'         indices for the sequence (`seq_ind`) and variant (`var_ind`).}
#'     \item{`set_names(names)`}{Set names for all variants.}
#'     \item{`rm_vars(var_names)`}{Remove one or more variants based on names in
#'         the `var_names` vector.}
#'
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
            private$check_var_ind(var_ind)
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
        extract_seq = function(var_ind, seq_ind) {
            private$check_ptr()
            private$check_var_ind(var_ind)
            private$check_seq_ind(seq_ind)
            return(view_var_genome_seq(self$genomes, var_ind - 1, seq_ind - 1))
        },


        # ----------*
        # __edit__ ----
        # ----------*

        # Change variant names
        set_names = function(names) {
            private$check_ptr()
            if (length(names) != self$n_vars()) {
                stop("names arg must be the same length as # variants", call. = FALSE)
            }
            var_inds <- 0:(length(names) - 1)
            set_var_set_var_names(self$genomes, var_inds, names)
            invisible(self)
        },

        # Remove one or more variants by name
        rm_vars = function(var_names) {
            private$check_ptr()
            self_names <- self$var_names()
            if (!all(var_names %in% self_names)) {
                stop("not all provided var_names are in this genome.", call. = FALSE)
            }
            if (anyDuplicated(var_names) != 0) {
                stop("one or more var_names are duplicate.", call. = FALSE)
            }
            var_inds <- match(var_names, self_names) - 1
            remove_var_set_vars(self$genomes, var_inds)
            invisible(self)
        },


        # Mutations:
        add_sub = function(var_ind, seq_ind, pos, nt) {
            private$check_pos(var_ind, seq_ind, pos)
            if (length(nt) != 1 | nchar(nt) != 1) {
                stop("nt arg must be a single character", call. = FALSE)
            }
            if (! nt %in% c("T", "C", "A", "G", "N")) {
                stop("nt arg must be one of \"T\", \"C\", \"A\", \"G\", or \"N\"",
                     call. = FALSE)
            }
            add_substitution(self$genomes, var_ind - 1, seq_ind - 1, nt, pos - 1)
            invisible(self)
        },

        add_ins = function(var_ind, seq_ind, pos, nts) {
            private$check_pos(var_ind, seq_ind, pos)
            if (length(nts) != 1) {
                stop("nts arg must be a single string", call. = FALSE)
            }
            if (! all(strsplit(nts, "")[[1]] %in% c("T", "C", "A", "G", "N"))) {
                stop("nts arg must only contain \"T\", \"C\", \"A\", \"G\", or \"N\"",
                     call. = FALSE)
            }
            add_insertion(self$genomes, var_ind - 1, seq_ind - 1, nts, pos - 1)
            invisible(self)
        },

        add_del = function(var_ind, seq_ind, pos, n_nts) {
            private$check_pos(var_ind, seq_ind, pos)
            if (!single_whole_number(n_nts)) {
                stop("n_nts arg must be a single whole number", call. = FALSE)
            }
            add_deletion(self$genomes, var_ind - 1, seq_ind - 1, n_nts, pos - 1)
            invisible(self)
        }

    ),


    private = list(
        # This should store a `XPtr<RefGenome>` to make sure it doesn't
        # go out of scope:
        reference = NULL,

        check_ptr = function() stopifnot(inherits(self$genomes, "externalptr")),

        check_seq_ind = function(seq_ind) {
            stopifnot(inherits(self$genomes, "externalptr"))
            if (seq_ind > self$n_seqs() | seq_ind < 1) {
                stop("seq_ind arg must be in range [1, <# sequences>]", call. = FALSE)
            }
        },
        check_var_ind = function(var_ind) {
            stopifnot(inherits(self$genomes, "externalptr"))
            if (var_ind > self$n_vars() | var_ind < 1) {
                stop("var_ind arg must be in range [1, <# variants>]", call. = FALSE)
            }
        },
        check_pos = function(var_ind, seq_ind, pos) {
            stopifnot(inherits(self$genomes, "externalptr"))
            if (seq_ind > self$n_seqs() | seq_ind < 1) {
                stop("seq_ind arg must be in range [1, <# sequences>]", call. = FALSE)
            }
            if (var_ind > self$n_vars() | var_ind < 1) {
                stop("var_ind arg must be in range [1, <# variants>]", call. = FALSE)
            }
            if (pos > self$sizes(var_ind)[seq_ind] | pos < 1) {
                stop("pos arg must be in range [1, <sequence size>]", call. = FALSE)
            }
        }
    )

)

variants$lock()





# >> region_generator class----
#' An R6 class representing information to generate probabilities by region.
#'
#'
#' @field distr_mean Mean of the distribution.
#' @field distr_sd Standard deviation of the distribution.
#' @field reg_len Length of each region.
#' @field n_reg Number of regions.
#' @field probs Vector of probabilities input manually.
#'
#' @return An object of class \code{region_generator}.
#'
#' @docType class
#'
#' @seealso \code{\link{reg_gen}}
#'
#' @export
#'
#' @format An \code{\link[R6]{R6Class}} generator object
#'
#' @describeIn reg_gen
#'
region_generator <- R6::R6Class(

    "region_generator",

    public = list(

        initialize = function(reg_len, n_reg, distr_mean, distr_sd, probs) {
            if (!single_whole_number(reg_len)) {
                stop("\nWhen initializing a `region_generator` object ",
                     "the `reg_len` field must be a single whole number.",
                     call. = FALSE)
            }
            if (!is.null(n_reg) && !single_whole_number(n_reg)) {
                stop("\nWhen initializing a `region_generator` object ",
                     "the `n_reg` field must be `NULL` or a single whole number.",
                     call. = FALSE)
            }
            if (!is.null(distr_mean) && !single_number(distr_mean)) {
                stop("\nWhen initializing a `region_generator` object ",
                     "the `distr_mean` field must be `NULL` or a single whole number.",
                     call. = FALSE)
            }
            if (!is.null(distr_sd) && !single_number(distr_sd)) {
                stop("\nWhen initializing a `region_generator` object ",
                     "the `distr_sd` field must be `NULL` or a single whole number.",
                     call. = FALSE)
            }
            if (!is.null(probs) && !inherits(probs, "numeric")) {
                stop("\nWhen initializing a `region_generator` object ",
                     "the `probs` field must be `NULL` or a numeric vector.",
                     call. = FALSE)
            }
            private$reg_len <- reg_len
            private$n_reg <- n_reg
            private$distr_mean <- distr_mean
            private$distr_sd <- distr_sd
            private$probs <- probs
        },

        print = function() {
            cat("< region_generator object >\n")
            for (x in c("reg_len", "n_reg", "distr_mean", "distr_sd", "probs")) {
                cat(sprintf("* %s field = %.3g\n", x, private[[x]]))
            }
            invisible(self)
        },

        generate = function(distr) {
            distr = match.arg(distr, c("gamma", "normal", "lognormal", "exponential",
                                       "uniform"))
            mu <- private$distr_mean
            sigma <- private$distr_sd
            n <- private$n_reg
            len <- private$reg_len
            if (is.null(n)) {
                stop("\nIf using a `region_generator` object to generate ",
                     "probabilities, you must provide it with the number ",
                     "of regions.",
                     call. = FALSE)
            }
            if (distr %in% c("gamma", "normal", "lognormal") &&
                (is.null(mu) || is.null(sigma))) {
                stop("\nIf using a `region_generator` object to generate ",
                     "probabilities from a gamma, normal, or lognormal distribution, ",
                     "you must provide it with a mean and SD.",
                     call. = FALSE)
            }
            if (distr == "exponential" && is.null(mu)) {
                stop("\nIf using a `region_generator` object to generate ",
                     "probabilities from an exponential distribution, ",
                     "you must provide it with a mean.",
                     call. = FALSE)
            }
            # Starting and ending points of each region:
            starts_ends <- cbind(sapply(0:(n-1), function(i) i * len),
                                 sapply(1:n, function(i) i * len))
            out_probs <- NULL
            if (distr == "gamma") {
                if (mu <= 0) {
                    stop("\nWhen using a gamma distribution, the mean must be > 0.",
                         call. = FALSE)
                }
                shape <- (mu / sigma)^2
                scale <- sigma^2 / mu
                out_probs <- apply(starts_ends, 1, function(x) {
                    pgamma(x[2], shape = shape, scale = scale) -
                        pgamma(x[1], shape = shape, scale = scale)
                })
                out_probs <- out_probs / sum(out_probs)
            } else if (distr == "normal") {
                out_probs <- apply(starts_ends, 1, function(x) {
                    pnorm(x[2], mu, sigma) - pnorm(x[1], mu, sigma)
                })
                out_probs <- out_probs / sum(out_probs)
            } else if (distr == "lognormal") {
                mulog <- log(mu / sqrt(1 + sigma / mu^2))
                sigmalog <- log(1 + sigma / mu^2)
                out_probs <- apply(starts_ends, 1, function(x) {
                    plnorm(x[2], meanlog = mulog, sdlog = sigmalog) -
                        plnorm(x[1], meanlog = mulog, sdlog = sigmalog)
                })
                out_probs <- out_probs / sum(out_probs)
            } else if (distr == "exponential") {
                rate_ <- 1 / mu
                out_probs <- apply(starts_ends, 1, function(x) {
                    pexp(x[2], rate = rate_) - pexp(x[1], rate = rate_)
                })
                out_probs <- out_probs / sum(out_probs)
            } else if (distr == "uniform") {
                out_probs <- rep(1, n) / n
            } else {
                stop(paste0("\nThe provided distribution (\"", distr,
                            "\") is not implemented."),
                     call. = FALSE)
            }

            return(out_probs)

        }
    ),

    private = list(

        reg_len = NULL,
        n_reg = NULL,
        distr_mean = NULL,
        distr_sd = NULL,
        probs = NULL

    )
)

region_generator$lock()


#' Create info to generate probabilities by region.
#'
#' @param reg_len Length of each region.
#' @param n_reg
#' @param distr_mean
#' @param distr_sd
#' @param probs
#'
#' @return
#' @export
#'
#' @examples
reg_probs <- function(reg_len = 1L,
                      n_reg = NULL,
                      distr_mean = NULL,
                      distr_sd = NULL,
                      probs = NULL) {
    region_generator$new(reg_len, n_reg, distr_mean, distr_sd, probs)
}
