//
// This file creates variants sequentially across scaffolds.
//

#include <RcppArmadillo.h>
#include <vector>  // begin, end
#include <algorithm>  // shuffle, lower_bound
#include <numeric>  // accumulate
#include <cmath>  // std::exp, std::log
#include <functional>  // std::function
#include <sitmo.h>  // sitmo::engine
#ifdef _OPENMP
#include <omp.h>  // omp
#endif



#include "gemino_types.h"  // integer types, SequenceSet, VariantSet
#include "vitter_algorithms.h"  // vitter_d


using namespace Rcpp;




// (this is faster than using std::pow(a, b))
// equivalent to a^b
inline double fast_pow(double a, double b) {
    return std::exp(b * std::log(a));
}

// Constants in use below
namespace variants {

    const std::string bases = "ACGT";

    // Cumulative probabilities of sampling an indel length from 1 to 10
    const std::vector<double> indel_probs_cumsum = {0.6321493, 0.864704, 0.9502561,
                                               0.9817289, 0.9933071, 0.9975665,
                                               0.9991335, 0.9997099, 0.999922,
                                               1};

    // Adding +1 to this so that when sampling, `engine() / variants::sitmo_max` is always
    // < 1. This is useful for discrete sampling bc I don't want there to be a number
    // that has a 1/2^32 chance of being sampled like what would happen if engine()
    // produced a number == `sitmo::prng_engine::max()`.
    const double sitmo_max = (double) sitmo::prng_engine::max() + 1.0;
}







// ======================================================================================
// ======================================================================================

//      Pre-iteration

// ======================================================================================
// ======================================================================================


// These functions are used for various computations before iterating through scaffolds.



// This function is to get the value of the two beta-distribution shapes that
// minimize the absolute difference between the expected and desired segregated-site
// divergence (i.e., mean pairwise differences)
// It is used in the R function `freq_probs` during the optimization step.
// [[Rcpp::export]]
double optim_prob(NumericVector v, NumericVector mean_pws_, NumericVector dens_,
                  double seg_div_) {
    double a = v[0];
    double b = v[1];
    NumericVector probs = dbeta(mean_pws_, a, b);
    probs = probs / dens_;
    double abs_diff = std::abs(seg_div_ - sum((probs / sum(probs)) * mean_pws_));
    return abs_diff;
}


//' Randomly choose scaffolds for segregating sites, weighted based on scaffold length.
//'
//' This function is used separately for indels and SNPs.
//'
//' The indices of the output matrix coincide with the order of scaffolds in the
//' \code{dna_set} input to \code{make_variants}.
//'
//' This function does NOT return an error if a scaffold is chosen more times
//' than its length.
//'
//' @param total_mutations The total number of mutations (SNPs and indels).
//' @param scaff_lens A vector of cumulative sums of scaffold lengths.
//'
//'
//' @return A numeric vector containing the number of mutations per scaffold.
//'
// [[Rcpp::export]]
std::vector<uint> sample_scaffs(const uint& total_mutations,
                         const std::vector<double>& scaff_lens_cumsum,
                         const std::vector<uint>& seeds) {

    const uint n_scaffs = scaff_lens_cumsum.size();
    const double sum_of_lens = scaff_lens_cumsum.back();
    const uint n_cores = seeds.size();
    // Creating output vector
    std::vector<uint> out_vec(n_scaffs, 0);

    #ifdef _OPENMP
    #pragma omp parallel default(shared) num_threads(n_cores) if(n_cores > 1)
    {
    #endif

    uint active_seed;

    // Write the active seed per core or just write one of the seeds.
    #ifdef _OPENMP
    uint active_thread = omp_get_thread_num();
    active_seed = seeds[active_thread];
    #else
    active_seed = seeds[0];
    #endif

    sitmo::prng_engine engine(active_seed);

    // Parallelize the Loop
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (unsigned int i = 0; i < total_mutations; i++){
        double rnd = ((double) engine() / variants::sitmo_max) * sum_of_lens;
        auto iter = lower_bound(scaff_lens_cumsum.begin(), scaff_lens_cumsum.end(), rnd);
        uint ind = distance(scaff_lens_cumsum.begin(), iter);
        #ifdef _OPENMP
        #pragma omp atomic
        #endif
        out_vec[ind] += 1;
    }

    #ifdef _OPENMP
    }
    #endif

    return out_vec;

}





// C++ equivalent of R's \code{choose} function.
//
// \emph{Note}: This function is not exported to R.
//
// @param n Integer value.
// @param k Integer value.
//
// @return Binomial coefficient (also integer).
//
// @name cpp_choose
//
inline unsigned cpp_choose(unsigned n, unsigned k) {
    if (k > n) {
        return 0;
    }
    if (k * 2 > n) {
        k = n - k;
    }
    if (k == 0) {
        return 1;
    }

    int result = n;
    for(int i = 2; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}




// Calculate mean pairwise differences between samples using a vector of nucleotide
//     frequencies.
//
//
// @param sample_segr Vector of nucleotide frequencies at a given segregating site for
//     all samples.
//
// @return Mean of the pairwise differences.
//
//
double cpp_mean_pairwise_freqs(const std::vector<int>& sample_segr) {

    int N_nt = sample_segr.size();

    if (N_nt != 4) {
        stop("Vector input to cpp_mean_pairwise_freqs must be of length 4.");
    }

    int N = accumulate(sample_segr.begin(), sample_segr.end(), 0.0);

    int total_pairs = cpp_choose(N, 2);
    double same_pairs = 0;
    for (int i = 0; i < N_nt; i++) {
        same_pairs += cpp_choose(sample_segr[i], 2);
    }

    double mean_pw = (total_pairs - same_pairs) / total_pairs;

    return mean_pw;

}



//' Get possible nucleotide distributions and their pairwise differences.
//'
//' Retrieve all combinations (with replacement) of nucleotide distributions that sum
//' to \code{N}, and, for each, calculate \eqn{\pi_{ji}}.
//'
//'
//' @param N Total number of individuals the frequencies must add to.
//'
//' @return List consisting of a matrix and a vector.
//'     The matrix (\code{List$combos}) contains all nucleotide frequencies that add
//'     to \code{N} (by row).
//'     The vector (\code{List$mean_pws}) contains the mean pairwise differences
//'     for a segregating site comprised of nucleotide frequencies present in each row
//'     of the matrix.
//'     For example, a segregating site for 10 haploid samples containing 3 As, 3 Cs,
//'     2 Gs, and 2 Ts would have a mean pairwise difference of 0.8222222.
//'
// [[Rcpp::export]]
List cpp_nt_freq(int N) {

    int a;
    std::vector<std::vector<int> > combos;
    std::vector<int> tmp_combos(4);
    std::vector<std::string> seq;
    std::vector<int> segr_bases(N);
    std::vector<double> mean_pws;
    double mean_pw;

    for (int d = 0; d <= floor(N / 4); d++) {
        for (int c = d; c <= floor((N - d) / 3); c++) {
            for (int b = c; b <= floor((N - c - d) / 2); b++) {
                a = N - b - c - d;
                tmp_combos[0] = a;
                tmp_combos[1] = b;
                tmp_combos[2] = c;
                tmp_combos[3] = d;
                combos.push_back(tmp_combos);
                mean_pw = cpp_mean_pairwise_freqs(tmp_combos);
                mean_pws.push_back(mean_pw);
                checkUserInterrupt();
            }
        }
    }

    arma::umat combo_mat(combos.size(), 4);

    for (int i = 0; i < combo_mat.n_rows; i++) {
        for (int j = 0; j < 4; j++) {
            combo_mat(i, j) = combos[i][j];
        }
    }

    return List::create(Named("combos") = combo_mat,
                        Named("mean_pws") = mean_pws);
}










// ======================================================================================
// ======================================================================================

//      Iteration

// ======================================================================================
// ======================================================================================


// These functions are used for iterating through scaffolds, ultimately creating a
// VariantSet object.






/*
 --------------------------

 InDels

 --------------------------
*/

// Sample indel lengths
uint indel_lengths(uint positions_left, sitmo::prng_engine& engine) {

    uint max_len = 10;
    if (positions_left < max_len) max_len = positions_left;
    const double sum_of_weight = variants::indel_probs_cumsum[(max_len - 1)];

    double rnd = ((double) engine() / variants::sitmo_max) * sum_of_weight;
    std::vector<double>::const_iterator iter = lower_bound(
        variants::indel_probs_cumsum.begin(),
        variants::indel_probs_cumsum.end(), rnd);
    // Here, `out` is from 0 --> (max_len - 1)
    uint out = distance(variants::indel_probs_cumsum.begin(), iter);
    out++; // Now `out` goes from 1-->max_len

    return out;
}


// Insertion nucleotides
inline std::vector<char> insertion_nucleos(const uint& length,
                                      const XPtr<SequenceSet>& reference,
                                      const uint& scaff_index,
                                      const uint& position,
                                      sitmo::prng_engine& engine) {

    std::vector<char> new_nucleos(1, reference->sequences[scaff_index][position]);
    // std::vector<char> new_nucleos(1, 'x');
    uint ind;
    for (uint i = 0; i < length; i++) {
        ind = ((double) engine() / variants::sitmo_max) * 4;
        new_nucleos.push_back(variants::bases[ind]);
    }

    return new_nucleos;
}

inline std::vector<char> deletion_nucleos(const uint& length) {
    std::vector<char> new_nucleos(length, '\0');
    return new_nucleos;
}



// Positions for insertions
inline std::vector<uint> insertion_sites(uint start_pos, uint length) {
    std::vector<uint> out(length + 1, start_pos);
    return out;
}


// Positions for deletions
inline std::vector<uint> deletion_sites(uint start_pos, uint length) {

    std::vector<uint> out(length);
    for (uint j = 0; j < length; j++) out[j] = start_pos + j;

    return out;
}





/*
 --------------------------

 SNPs

 --------------------------
*/


// Make snp output char vector
inline std::vector<char> snp_nucleos(const arma::umat& snp_combo_mat,
                                const std::vector<double>& snp_probs_cumsum,
                                sitmo::prng_engine& engine) {

    const double sum_of_weight = snp_probs_cumsum.back();

    double rnd = ((double) engine() / variants::sitmo_max) * sum_of_weight;
    std::vector<double>::const_iterator iter = lower_bound(snp_probs_cumsum.begin(),
                                                      snp_probs_cumsum.end(), rnd);
    uint ind = distance(snp_probs_cumsum.begin(), iter);
    std::vector<uint> combo = arma::conv_to<std::vector<uint>>::from(snp_combo_mat.row(ind));
    shuffle(combo.begin(), combo.end(), engine);

    std::vector<char> seq_out;

    for (uint i = 0; i < 4; i++) {
        for (uint j = 0; j < combo[i]; j++) {
            seq_out.push_back(variants::bases[i]);
        }
    }
    shuffle(seq_out.begin(), seq_out.end(), engine);

    return seq_out;

}






/*
 --------------------------

 Putting it together

 --------------------------
*/

// Change variant_set and return length for one segregating site
uint one_site(
        XPtr <VariantSet>& variant_set,
        const sint& n, const uint& N, const uint& S, const uint& current_pos,
        const uint& n_vars, const uint& scaff_len, const uint& scaff_index,
        const XPtr<SequenceSet>& reference,
        sitmo::prng_engine& engine,
        const arma::umat& snp_combo_mat,
        const std::vector<double>& snp_probs_cumsum,
        const double& snp_p,
        const double& insertion_p
    ) {

    uint length;

    // Number of variants to get the indel, and random indices of that # of variants
    uint n_vars_w_ins;
    std::vector<uint> vars_w_ins;

    // Storing nucleos, sites, and scaffold modifier
    std::vector<char> nucleos;
    std::vector<uint> sites;
    sint scaff_mod;

    // Random number (0,1) indicating which type of segregating site (SNP, insertion,
    // or deletion) a site will be
    double type_rnd = (double) engine() / variants::sitmo_max;

    /* ~~~ SNP ~~~ */
    if (type_rnd < snp_p) {
        length = 1;
        nucleos = snp_nucleos(snp_combo_mat, snp_probs_cumsum, engine);
        for (uint v = 0; v < n_vars; v++) {
            variant_set->variant_info[v].nucleos[scaff_index].push_back(nucleos[v]);
            variant_set->variant_info[v].sites[scaff_index].push_back(current_pos);
        }
    /* ~~~ InDel ~~~ */
    } else {
        /* --- Insertion --- */
        if (type_rnd < (snp_p + (1 - snp_p) * insertion_p)) {
            // Insertions do not need to be length-limited
            length = indel_lengths(999, engine);
            nucleos = insertion_nucleos(length, reference, scaff_index,
                                        current_pos, engine);
            sites = insertion_sites(current_pos, length);
            scaff_mod = length;
            // Changing back to 1 bc insertions don't need to be skipped over like
            // deletions bc they don't take up existing nucleotides
            length = 1;
        /* --- Deletion --- */
        } else {
            // Deletions need to be length-limited, else we may run out of scaffold
            length = indel_lengths(N - S - n + 1, engine);
            nucleos = deletion_nucleos(length);
            sites = deletion_sites(current_pos, length);
            scaff_mod = -1 * length;
        }
        n_vars_w_ins = ((double) engine() / variants::sitmo_max) * (n_vars - 1);
        n_vars_w_ins++;
        vars_w_ins = vitter_d(n_vars_w_ins, n_vars, engine);
        for (uint v : vars_w_ins) {
            for (uint i = 0; i < nucleos.size(); i++) {
                variant_set->variant_info[v].nucleos[scaff_index].push_back(nucleos[i]);
                variant_set->variant_info[v].sites[scaff_index].push_back(sites[i]);
            }
            variant_set->variant_info[v].scaffold_lengths[scaff_index] += scaff_mod;
        }
    }
    return length;
}



// Iterate through one scaffold, changing segregating each segregating site with
// `one_site`

void one_scaff(
        XPtr <VariantSet>& variant_set,
        const uint& n_segr,
        const uint& scaff_len,
        const uint& scaff_index,
        const XPtr<SequenceSet>& reference,
        sitmo::prng_engine& engine,
        const arma::umat& snp_combo_mat,
        const std::vector<double>& snp_probs_cumsum,
        double snp_p = 0.9,
        double insertion_p = 0.5,
        double n2N = 50,
        double alpha = 0.8
    ) {

    uint n_vars = arma::sum(snp_combo_mat.row(0));

    // These values are copied bc n and N will be changing
    sint n = n_segr;
    uint N = scaff_len;

    // Commented this out bc this will crash R if run in parallel and stop happens.
    // if (alpha > 1 || alpha < 0) stop("Invalid alpha. It must be (0,1).");

    // The # positions to skip before taking the next one (0 to (N - n - 1))
    uint S;
    // Keeping track of the current position
    sint64 current_pos = -1;  // (starts at -1 so it can reach 0 on the first skip)

    uint length; // Length of segtregating sites (SNPs are obviously always 1)

    // Stores function to sequentially find `S`
    std::function<uint(const sint&, const uint&, sitmo::prng_engine&,
                       const double)> algorithm;

    if ((fast_pow(n, 2) / N) > n2N) {
        algorithm = algorithm_d2_S;
    } else {
        algorithm = algorithm_d1_S;
    }
    while (n > 0) {
        if (n > 1) {
            S = algorithm(n, N, engine, alpha);
        } else { // At n = 1, D2 divides by zero, but below works just fine
            S = ((double) engine() / variants::sitmo_max) * N;
        }
        current_pos += S + 1;
        // This function returns the length, but modifies variant_set
        length = one_site(variant_set, n, N, S, current_pos, n_vars, scaff_len,
                          scaff_index, reference, engine, snp_combo_mat, snp_probs_cumsum,
                          snp_p, insertion_p);
        current_pos += (length - 1);
        n--;
        N -= (S + length);
    }

    return;
}





//' Inner function to create a C++ \code{VariantSet} object
//'
//' A \code{VariantSet} object constitutes the majority of information in a
//' \code{variants} object (other than the reference genome) and is located in
//' the \code{variant_set} field.
//'
//' @param n_mutations Integer vector of the total number of mutations (SNPs or indels)
//'     for each scaffold.
//' @param reference External pointer to a C++ \code{SequenceSet} object that
//'     represents the reference genome.
//' @param snp_combo_mat Matrix of all possible nucleotide combinations among all
//'     variants per SNP.
//' @param snp_probs_cumsum Vector of sampling probabilities for each row in
//'     \code{snp_combo_mat}.
//' @param seeds Vector of seeds, the length of which dictates how many cores will be
//'     used.
//' @param snp_p Proportion of mutations that are SNPs. Defaults to 0.9.
//' @param insertion_p Proportion of \emph{indels} that are insertions. Defaults to 0.5.
//' @param n2N A numeric threshold placed on the algorithm used to find new locations.
//'     This is not recommended to be changed. Defaults to 50.
//' @param alpha A numeric threshold placed on the algorithm used to find new locations.
//'     This is not recommended to be changed. Defaults to 0.8.
//'
//'
//' @return An external pointer to a \code{VariantSet} object in C++.
//'
//[[Rcpp::export]]
XPtr<VariantSet> make_variant_set(
        const std::vector<uint>& n_mutations,
        const XPtr<SequenceSet>& reference,
        const arma::umat& snp_combo_mat,
        const std::vector<double>& snp_probs_cumsum,
        std::vector<uint> seeds,
        double snp_p = 0.9,
        double insertion_p = 0.5,
        double n2N = 50,
        double alpha = 0.8
    ) {

    const std::vector<uint> scaff_lens = reference->seq_sizes;
    const uint n_scaffs = scaff_lens.size();
    const uint n_cores = seeds.size();
    const uint n_vars = arma::sum(snp_combo_mat.row(0));

    if (n_mutations.size() != n_scaffs) stop("n_mutations is incorrect length.");
    if (alpha > 1 || alpha < 0) stop("Invalid alpha. It must be [0,1].");
    if (snp_p > 1 || snp_p < 0) stop("Invalid snp_p. It must be [0,1].");
    if (insertion_p > 1 || insertion_p < 0) {
        stop("Invalid insertion_p. It must be [0,1].");
    }

    uint total_segrs = accumulate(n_mutations.begin(), n_mutations.end(), 0.0);

    XPtr<VariantSet> variant_set(new VariantSet(n_vars, total_segrs, scaff_lens), true);

    #ifdef _OPENMP
    #pragma omp parallel shared(variant_set) num_threads(n_cores) if (n_cores > 1)
    {
    #endif

    uint active_seed;

    // Write the active seed per core or just write one of the seeds.
    #ifdef _OPENMP
    active_seed = seeds[omp_get_thread_num()];
    #else
    active_seed = seeds[0];
    #endif

    sitmo::prng_engine engine(active_seed);

    // Parallelize the Loop
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint s = 0; s < n_scaffs; s++) {
        one_scaff(variant_set, n_mutations[s], scaff_lens[s], s, reference,
                  engine, snp_combo_mat, snp_probs_cumsum,
                  snp_p, insertion_p, n2N, alpha);
    }
    #ifdef _OPENMP
    }
    #endif

    return variant_set;
}


