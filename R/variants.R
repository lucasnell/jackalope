
#' Compute nucleotide-frequency probabilities.
#'
#' This function computes the sampling probabilities necessary to conduct weighted
#' sampling from the vector of mean pairwise differences and get an overall mean
#' divergence for all segregating sites equal to one that is input to the
#' function.
#'
#'
#' @param mean_pws A vector of mean pairwise differences for all potential nucleotide
#'     frequencies.
#' @param seg_div The overall segregating site divergence that is desired.
#' @param threshold Threshold for how close you want the output divergence to be in
#'     relation to the desired one.
#' @param make_converge Boolean for whether to return an error if it does not converge
#'     successfully. Defaults to \code{TRUE}.
#'
#' @return A vector of probabilities of the same length as \code{mean_pws}.
#'
#'
freq_probs <- function(mean_pws, seg_div, threshold = 0.0001, make_converge = TRUE) {
    # Kernel density estimation
    dens_obj <- density(mean_pws)
    # Assigning densities to each value in mean_pws
    dens <- sapply(mean_pws, function(x) {
        d <- dens_obj$y[abs(x - dens_obj$x) == min(abs(x - dens_obj$x))]
        if (length(d) > 1) return(mean(d))
        return(d)
    })
    # This returns the value of the two beta-distribution shapes that
    # minimize the absolute difference between the expected and desired segregated-site
    # divergence (i.e., mean pairwise differences)
    # It also returns the absolute difference at those values.
    # The function optim_prob is a C++ function for minimization.
    optim_out <- suppressWarnings(optim(par = list(1, 1), fn = optim_prob,
                                        mean_pws_ = mean_pws, dens_ = dens,
                                        seg_div_ = seg_div, method = "Nelder-Mead"))

    # Checks on optim convergence
    if (make_converge & optim_out$convergence == 1) {
        stop("The iteration limit maxit had been reached.")
    }
    if (make_converge & optim_out$convergence == 10) {
        stop("Degeneracy of the Nelder–Mead simplex.")
    }
    if (optim_out$value > threshold) {
        stop('Optimization did not find a sufficiently accurate value.')
    }

    probs <- dbeta(mean_pws, optim_out$par[1], optim_out$par[2]) / dens

    probs <- probs / sum(probs)

    return(probs)

}



#' Return SNP nucleotide combinations and their sampling weights.
#'
#'
#' @param n_vars The number of variants.
#' @param seg_div The desired segregating-site divergence.
#' @param snp_site_prop The proportion of segregating sites representing SNPs (vs indels).
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
get_snp_combos_weights <- function(n_vars, seg_div, snp_site_prop) {

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
    nt_probs <- freq_probs(nt_freq$mean_pws, snp_div);
    # My weighted-sampling algorithm uses cumulative sums of weights
    nt_probs <- cumsum(nt_probs)

    return(list(combo_mat = nt_freq$combos, probs_cumsum = nt_probs));
}









#' Construct a \code{variants} object from a reference genome and summary statistics.
#'
#' This function creates a \code{\link{variants}} object, which is designed to be a
#' low-memory way to store variants.
#' (The variants class also prevents what might otherwise be an annoyingly long list
#' from ever printing in the console.)
#'
#' @param dna_set_in A \code{dna_set} object of sequences representing the reference
#'     genome.
#' @param n_vars The number of variants to create.
#' @param seg_prop The proportion of sites in the genome that are segregating. Defaults
#'     to 0.01414.
#' @param seg_div The mean pairwise divergence at segregating sites. Defaults to 0.7072.
#' @param n_cores Number of cores to use. Defaults to 1.
#' @param snp_prop The proportion of mutations (not sites) that are SNPs. Defaults to 0.9.
#' @param insertion_prop The proportion of mutations (not sites) that are indels.
#'     Defaults to 0.5.
#' @param n2N A numeric threshold placed on the algorithm used to find new locations.
#'     This is not recommended to be changed. Defaults to 50.
#' @param alpha A numeric threshold placed on the algorithm used to find new locations.
#'     This is not recommended to be changed. Defaults to 0.8.
#'
#' @return A \code{\link{variants}} object.
#'
#' @seealso \code{\link{variants}}
#'
#' @export
#'
#'
#' @examples
#' n_vars <- 10
#' dna_set_in <- dna_set$new(rando_seqs(100, 100))
#' set.seed(1)
#' varseq_out <- make_variants(dna_set_in, n_vars)
#'
make_variants <- function(dna_set_in, n_vars, theta_w = 0.0050, theta_pi = 0.0045,
                          n_cores = 1, snp_prop = 0.9, insertion_prop = 0.5,
                          n2N = 50, alpha = 0.8) {

    # Proportion of sites that are segregating
    seg_prop = theta_w * sum(1 / 1:(n_vars-1))
    # Mean divergence at segregating sites
    seg_div <- { theta_pi * n_vars^2 } / { choose(n_vars, 2) * seg_prop }


    # Useful info from the input dna_set
    scaff_lens <- seq_sizes_SequenceSet(dna_set_in$sequence_set)
    total_seg <- round(sum(scaff_lens) * seg_prop)
    n_scaffs <- length(scaff_lens)

    # Proportions of segregating sites represented by SNPs and indels
    # (Accounts for the fact that the average indel length (1.581523) is > 1bp)
    snp_site_prop <- snp_prop / ((1 - snp_prop) * 1.581523 + snp_prop)

    # Total number of mutations (SNPs and indels; not on a site basis)
    # for all scaffolds
    total_mutations <- round((snp_site_prop / snp_prop) * total_seg)

    # SNP nucleotide combinations and their sampling weights.
    snp_combos_weights <- get_snp_combos_weights(n_vars, seg_div, snp_site_prop)

    # Setting seeds for thread-safe C++ pseudo-random number generators (1 per core)
    seeds <- sample.int(2^31 - 1, n_cores)

    # Sampling the number of mutations per scaffold
    n_mutations <- sample_scaffs(total_mutations, cumsum(scaff_lens), seeds)

    # Setting new seeds for the next step
    seeds <- sample.int(2^31 - 1, n_cores)

    variant_set <- make_variant_set(
        n_mutations,
        dna_set_in$sequence_set,
        snp_combos_weights$combo_mat, snp_combos_weights$probs_cumsum,
        seeds, snp_prop, insertion_prop, n2N, alpha)


    var_obj <- variants$new(dna_set_in$sequence_set, variant_set)

    return(var_obj)
}




# Below study simulates 1:9 ratio of indels to SNPs:
# Albers, C. A., G. Lunter, D. G. MacArthur, G. McVean, W. H. Ouwehand, and R. Durbin.
#     2011. Dindel: accurate indel calls from short-read data. Genome Research 21:961–973.
# Also...
#   The length of the indels varied from 1 to 10 bp, and the length distribution was
#   such that the number of indels with length l was proportional to exp(-l), similar to
#   what is observed in real data sets.
