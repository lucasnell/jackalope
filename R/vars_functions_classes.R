

# vars_ssites_info ----
#' An R6 class representing information for ssites method.
#'
#' @noRd
#'
#' @importFrom R6 R6Class
#'
vars_ssites_info <- R6Class(

    "vars_ssites_info",

    public = list(

        initialize = function(mats) {

            extra_msg <- paste(" Please only create these objects using the vars_ssites",
                               "function, NOT using vars_ssites_info$new().")
            if (!inherits(mats, "list") || !all(sapply(mats, is.numeric)) ||
                !all(sapply(mats, inherits, what = "matrix")) ||
                any(sapply(mats, function(x) any(x < 0)))) {
                stop("\nWhen initializing a vars_ssites_info object, you need to use ",
                     "a list of numeric matrices of values >= 0.", extra_msg,
                     call. = FALSE)
            }

            private$r_mats <- mats
        },

        print = function(...) {

            cat("< Seg. site variant-creation info >\n")
            cat(sprintf("# Number of variants = %i\n", ncol(private$r_mats[[1]]) - 1))
            cat(sprintf("# Number of sites = %s\n", format(as.integer(sum(sapply(
                private$r_mats, nrow))), big.mark = ",")))
            invisible(self)

        },

        mats = function() return(private$r_mats)

    ),

    private = list(

        r_mats = NULL

    ),

    lock_class = TRUE

)




# vars_vcf_info ----
#' An R6 class representing information for VCF method.
#'
#' @noRd
#'
#' @importFrom R6 R6Class
#'
vars_vcf_info <- R6Class(

    "vars_vcf_info",

    public = list(

        initialize = function(fn,
                              print_names) {

            err <- function(a, b) {
                msg <- paste0("\nWhen initializing a vars_vcf_info object, the ",
                              "argument `", a, "` should be ", b, ". ",
                              "Please only create these objects using the vars_vcf ",
                              "function, NOT using vars_vcf_info$new().")
                stop(msg, call. = FALSE)
            }

            if (!is_type(fn, "character", 1L)) {
                err("fn", "a single character")
            }
            if (!is_type(print_names, "logical", 1L)) {
                err("print_names", "a single logical")
            }

            private$r_fn <- fn
            private$r_print_names <- print_names
        },

        print = function(...) {

            digits <- max(3, getOption("digits") - 3)

            cat("< VCF variant-creation info >\n")
            cat(sprintf("# File name: %s\n", fn))

            invisible(self)

        },

        fn = function() return(private$r_fn),
        print_names = function() return(private$r_print_names)

    ),

    private = list(

        r_fn = NULL,
        r_print_names = NULL

    ),

    lock_class = TRUE

)




# vars_phylo_info ----
#' An R6 class representing information for phylo method.
#'
#' @noRd
#'
#' @importFrom R6 R6Class
#'
vars_phylo_info <- R6Class(

    "vars_phylo_info",

    public = list(

        initialize = function(phylo) {

            extra_msg <- paste(" Please only create these objects using the vars_phylo",
                               "function, NOT using vars_phylo_info$new().")
            if (!inherits(phylo, "list") ||
                !all(sapply(phylo, inherits, what = "phylo"))) {
                stop("\nWhen initializing a vars_phylo_info object, you need to use",
                     "a list of `phylo` object(s).", extra_msg, call. = FALSE)
            }

            private$r_phylo <- phylo
        },

        print = function(...) {

            cat("< Phylo variant-creation info >\n")
            cat(sprintf("# Number of variants = %i\n",
                        length(private$r_phylo[[1]]$tip.label)))
            cat(sprintf("# Number of trees = %i\n", length(private$r_phylo)))

            invisible(self)

        },

        phylo = function() return(private$r_phylo)

    ),

    private = list(

        r_phylo = NULL

    ),

    lock_class = TRUE

)






# vars_theta_info ----
#' An R6 class representing information for theta method.
#'
#' @noRd
#'
#' @importFrom R6 R6Class
#'
vars_theta_info <- R6Class(

    "vars_theta_info",

    public = list(

        initialize = function(phylo, theta) {

            extra_msg <- paste(" Please only create these objects using the vars_theta",
                               "function, NOT using vars_theta_info$new().")
            if (!inherits(phylo, "phylo") || !single_number(theta) || theta < 0) {
                stop("\nWhen initializing a vars_theta_info object, you need to use",
                     "a `phylo` object and a single number >= 0.", extra_msg,
                     call. = FALSE)
            }

            private$r_phylo <- phylo
            private$r_theta <- theta
        },

        print = function(...) {

            digits <- max(3, getOption("digits") - 3)

            cat("< Theta variant-creation info >\n")
            cat(sprintf(sprintf("# Theta: %%.%ig\n", digits), private$r_theta))
            cat("# Phylogenetic tree:\n")
            print(private$r_phylo)

            invisible(self)

        },

        phylo = function() return(private$r_phylo),
        theta = function() return(private$r_theta)

    ),

    private = list(

        r_phylo = NULL,
        r_theta = NULL

    ),

    lock_class = TRUE

)




# vars_gtrees_info ----
#' An R6 class representing information for gtrees method.
#'
#' @noRd
#'
#' @importFrom R6 R6Class
#'
vars_gtrees_info <- R6Class(

    "vars_gtrees_info",

    public = list(

        initialize = function(trees) {

            extra_msg <- paste(" Please only create these objects using the vars_gtrees",
                               "function, NOT using vars_gtrees_info$new().")
            if (!inherits(trees, "list") ||
                !all(sapply(trees, inherits, what = "character"))) {
                stop("\nWhen initializing a vars_gtrees_info object, you need to use",
                     "a list of character vectors.", extra_msg, call. = FALSE)
            }

            private$r_trees <- trees
        },

        print = function(...) {

            digits <- max(3, getOption("digits") - 3)

            cat("< Gene trees variant-creation info >\n")
            cat(sprintf("# Number of chromosomes: %i\n", length(private$r_trees)))
            cat(sprintf("# Number of variants: %i\n",
                        length(ape::read.tree(text = private$r_trees[[1]][1])$tip.label)))
            cat(sprintf("# Total trees: %i\n", sum(sapply(private$r_trees, length))))

            invisible(self)

        },

        trees = function() return(private$r_trees)

    ),

    private = list(

        r_trees = NULL

    ),

    lock_class = TRUE

)


