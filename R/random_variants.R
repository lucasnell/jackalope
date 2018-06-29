

#' Return SNP nucleotide combinations and their sampling weights.
#'
#'
#' @param n_vars The number of variants.
#' @param seg_div The desired segregating-site divergence.
#' @param snp_site_prop The proportion of segregating sites representing SNPs (vs indels).
#' @param threshold Threshold for how close you want the output divergence to be in
#'     relation to the desired one.
#' @param make_converge Boolean for whether to return an error if it does not converge
#'     successfully. Defaults to \code{TRUE}.
#' @param optim_opts List of additional arguments to pass to `stats::optim`.
#'
#'
#' @return A list containing a matrix and numeric vector.
#' \describe{
#'     \item{matrix}{The \code{combo_mat} field contains all combinations of nucleotides
#'     that sum to \code{n_vars}.
#'     Each row coincides with a single SNP, and each cell is the frequency of
#'     a given nucleotide (A, C, G, or T) among all variants at that SNP; the order
#'     of columns in relation to specific nucleotides is arbitrary.}
#'     \item{numeric}{The \code{probs_cumsum} contains the cumulative sum of sampling
#'     probabilities for each row in the \code{combo_mat} field matrix.}
#' }
#'
#'
#' @noRd
#'
#'
get_snp_combos_weights <- function(n_vars, seg_div, snp_site_prop,
                                   threshold = 0.0001, make_converge = TRUE,
                                   optim_opts = list()) {

    indel_site_prop <- 1 - snp_site_prop

    # List object containing...
    #     1. All possible nucleotide distributions for n_vars variants
    #     2. Pairwise differences for all of these combinations
    nt_freq <- cpp_nt_freq(n_vars)
    # Because indels only have two states, their pairwise divergence maxes out
    # a little over 0.5 (depending on n_vars)
    all_pairs <- choose(n_vars, 2)
    indel_pw <- sapply(1:(n_vars-1),
                       function(v) {
                           same_pairs <- choose(v, 2) + choose(n_vars - v, 2)
                           return((all_pairs - same_pairs) / all_pairs)
                       })
    # I'm uniformly choosing from 1 to (n_vars - 1) variants to get each indel, so the
    # mean pairwise divergence for that distribution is probably less than seg_div.
    # I need to adjust the snp divergence accordingly
    mean_indel_pw <- mean(indel_pw)
    snp_div <- (seg_div - mean_indel_pw * indel_site_prop) / snp_site_prop



    # Compute probabilities of being sampled for each pairwise difference.
    # These probs were designed to, on average, get the desired SNP
    # segregating site divergence (`snp_div`).
    # nt_probs <- freq_probs(nt_freq$mean_pws, snp_div, threshold, make_converge,
    #                        optim_opts)
    # Kernel density estimation
    dens_obj <- stats::density(nt_freq$mean_pws)
    # Assigning densities to each value in nt_freq$mean_pws
    dens <- sapply(nt_freq$mean_pws, function(x) {
        d <- dens_obj$y[abs(x - dens_obj$x) == min(abs(x - dens_obj$x))]
        if (length(d) > 1) return(mean(d))
        return(d)
    })
    # This returns the value of the two beta-distribution shapes that
    # minimize the absolute difference between the expected and desired segregated-site
    # divergence (i.e., mean pairwise differences)
    # It also returns the absolute difference at those values.
    # The function optim_prob is a C++ function for minimization.
    args <- c(optim_opts, list(par = list(1, 1), fn = quote(optim_prob),
                               mean_pws_ = quote(nt_freq$mean_pws), dens_ = quote(dens),
                               seg_div_ = seg_div))
    args[['method']] <- "Nelder-Mead"
    optim_out <- suppressWarnings(do.call(stats::optim, args))

    # Checks on stats::optim convergence
    if (make_converge & optim_out$convergence == 1) {
        stop("The iteration limit maxit had been reached.")
    }
    if (make_converge & optim_out$convergence == 10) {
        stop("Degeneracy of the Nelder-Mead simplex.")
    }
    if (optim_out$value > threshold) {
        stop('Optimization did not find a sufficiently accurate value.')
    }

    probs <- stats::dbeta(nt_freq$mean_pws, optim_out$par[1], optim_out$par[2]) / dens

    nt_probs <- probs / sum(probs)

    return(list(combo_mat = nt_freq$combos, probs = nt_probs));
}









#' Construct a \code{VarSet} object from a reference genome and summary statistics.
#'
#' This function creates a \code{VarSet} object, which is designed to be a
#' low-memory way to store variants.
#' (The variants class also prevents what might otherwise be an annoyingly long list
#' from ever printing in the console.)
#'
#' @param dna_set_in A \code{dna_set} object of sequences representing the reference
#'     genome.
#' @param n_vars The number of variants to create.
#' @param theta_w Watterson's estimator for the focal population.
#' @param theta_pi Average nucleotide diversity for the focal population.
#' @param snp_probs Relative probabilities of substitution types:
#'     "A", "C", "G", and "T" respectively. Defaults to `rep(0.25, 4)`.
#' @param indel_probs Relative probabilities of indel types and sizes.
#'     If insertions and deletions have the same probabilities, this is simply
#'     a numeric vector where the value in location `i` indicates the relative
#'     probability of an indel of size `i`.
#'     If insertions and deletions do not have the same probabilities, then
#'     this argument should be a list where the `insertion` and `deletion` fields are
#'     numeric vectors specifying relative probabilities for insertions and deletions,
#'     respectively.
#'     Note that if specifying a list, the proportion of insertions to deletions will
#'     be `sum(indel_probs$insertions) / sum(indel_probs$deletions)`.
#'     Defaults to `exp(-1:-10)`.
#' @param snp_proportion The proportion of mutations (not sites) that are SNPs.
#'     Defaults to the proportion calculated using the average ratio of indels
#'     to substitutions in eukaryotes from Sung et al. (2016).
#' @param n_cores Number of cores to use. Defaults to 1.
#' @param n2N A numeric threshold placed on the algorithm used to find new locations.
#'     This is not recommended to be changed. Defaults to 50.
#' @param alpha A numeric threshold placed on the algorithm used to find new locations.
#'     This is not recommended to be changed. Defaults to 0.8.
#'
#' @return A \code{VarSet} object.
#'
#'
#' @references
#' Sung, W., M. S. Ackerman, M. M. Dillon, T. G. Platt, C. Fuqua, V. S. Cooper, and M.
#' Lynch. 2016.
#' Evolution of the insertion-deletion mutation rate across the tree of life.
#' *G3: Genes|Genomes|Genetics* __6__:2583-2591.
#'
#' @examples
#' \dontrun{
#' n_vars <- 10
#' dna_set_in <- dna_set$new(gemino:::rando_seqs(100, 100))
#' set.seed(1)
#' varseq_out <- random_variants(dna_set_in, n_vars,
#'                               theta_w = 0.0045, theta_pi = 0.005)
#' }
#'
random_variants <- function(dna_set_in, n_vars, theta_w, theta_pi,
                            indel_probs = exp(-1:-10),
                            snp_probs = rep(0.25, 4),
                            snp_proportion = NULL,
                            n_cores = 1, n2N = 50, alpha = 0.8) {

    # If nothing provided, use data from Sung et al. (2016):
    if (is.null(snp_proportion)) {
        rates_ <- gemino::evo_rates
        snp_proportion <- rates_$subs[rates_$domain == 'Eukarya'] /
            rates_$indels[rates_$domain == 'Eukarya']
        snp_proportion <- mean(snp_proportion)
        snp_proportion <- snp_proportion / (1 + snp_proportion)
    }

    # Proportion of sites that are segregating
    seg_prop = theta_w * sum(1 / 1:(n_vars-1))
    # Mean divergence at segregating sites
    seg_div <- { theta_pi * n_vars^2 } / { choose(n_vars, 2) * seg_prop }


    # Useful info from the input dna_set
    n_seqs <- see_ref_n_seq(dna_set_in$sequence_set)
    seq_lens <- sapply(0:(n_seqs-1), see_ref_seq_size, ref_ = dna_set_in$sequence_set)
    total_seg <- round(sum(seq_lens) * seg_prop)

    # Standardize and get info from `indel_props`
    if (inherits(indel_probs, "numeric")) {
        indel_probs <- indel_probs / (2 * sum(indel_probs))
        indel_probs <- list(insertions = indel_probs, deletions = indel_probs)
    } else if (inherits(indel_probs, "list")) {
        if (!all(names(indel_probs) %in% c("insertions", "deletions")) |
            length(indel_probs) != 2) {
            stop("\nindel_probs argument to random_variants function must be a list ",
                 "of length two with names \"insertions\", \"deletions\".",
                 call. = FALSE)
        }
        ip_sum <- sum(indel_probs$insertions + indel_probs$deletions)
        indel_probs$insertions <- indel_probs$insertions / ip_sum
        indel_probs$deletions <- indel_probs$deletions / ip_sum
    } else {
        stop("\nindel_probs argument to random_variants function must be a ",
             "numeric vector or list", call. = FALSE)
    }
    if (max(sapply(indel_probs, length)) >= 100) {
        warning("`random_variants` is not designed for large indels. ",
                "I recommend another of the variant-creation methods.",
                call. = FALSE)
    }
    # Average indel size:
    avg_indel <- sum(sapply(indel_probs, function(x) 1:length(x) * x))
    # Now make these into the probabilities that will be used alongside `snp_proportion`
    indel_probs <- lapply(indel_probs, function(x) x * (1 - snp_proportion))

    # Standardize `snp_props`
    if (length(snp_probs) != 4) {
        stop("\nsnp_probs argument to random_variants must be of length 4", call. = FALSE)
    }
    snp_probs <- snp_probs / sum(snp_probs)

    # Average mutation size
    # (Accounts for the fact that the average indel length `avg_indel` is > 1bp)
    avg_mutation <- avg_indel * (1 - snp_proportion) + snp_proportion

    # Total number of mutations (SNPs and indels; not on a site basis)
    # for all sequences
    total_mutations <- round(total_seg / avg_mutation)

    # Proportions of segregating sites represented by SNPs (versus indels)
    snp_site_prop <- snp_proportion / ((1 - snp_proportion) * avg_indel + snp_proportion)

    # SNP nucleotide combinations and their sampling weights.
    snp_combos_weights <- get_snp_combos_weights(n_vars, seg_div, snp_site_prop,
                                                 threshold = 40)

    snp_combo_list <- split(snp_combos_weights$combo_mat, row(snp_combos_weights$combo_mat))

    # Column 1 is sampling weight, column 2 is type of mutation
    sampling_weights <- c(
        list(cbind(snp_combos_weights$probs * snp_proportion, 0)),
        lapply(c("insertions", "deletions"), function(x) {
            cbind(indel_probs[[x]], ifelse(x == "insertions", 1, 2))
        }))
    sampling_weights <- do.call(rbind, sampling_weights)

    mutation_sizes <- c(rep(0, sum(sampling_weights[,2] == 0)),
                        1:sum(sampling_weights[,2] == 1),
                        1:sum(sampling_weights[,2] == 2))

    # Sampling the number of mutations per sequence
    n_mutations <- sample_seqs(total_mutations, seq_lens, n_cores)


    var_set <- make_variants_(n_mutations, dna_set_in$sequence_set, snp_combo_list,
                              mutation_probs = sampling_weights[,1],
                              mutation_types = sampling_weights[,2],
                              mutation_sizes = mutation_sizes, n_cores = n_cores,
                              n2N = 50, alpha = 0.8)

    # var_obj <- variants$new(dna_set_in$sequence_set, variant_set)
    var_obj <- var_set

    return(var_obj)
}




# Below study simulates 1:9 ratio of indels to SNPs:
# Albers, C. A., G. Lunter, D. G. MacArthur, G. McVean, W. H. Ouwehand, and R. Durbin.
#     2011. Dindel: accurate indel calls from short-read data. Genome Research 21:961-973.
# Also...
#   The length of the indels varied from 1 to 10 bp, and the length distribution was
#   such that the number of indels with length l was proportional to exp(-l), similar to
#   what is observed in real data sets.
