/*
 ********************************************************

 This file creates variants randomly across sequences.

 ********************************************************
 */

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



#include "gemino_types.h"  // integer types
#include "sequence_classes.h"  // Ref* and Var* classes
#include "vitter_algorithms.h"  // vitter_d
#include "alias.h" // alias sampling


using namespace Rcpp;





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






/*
 ======================================================================================
 ======================================================================================

 Pre-iteration
 These functions are used for various computations before iterating through sequences.

 ======================================================================================
 ======================================================================================
 */



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


//' Randomly choose sequences for segregating sites, weighted based on sequence length.
//'
//' This function is used separately for indels and SNPs.
//'
//' The indices of the output matrix coincide with the order of sequences in the
//' \code{dna_set} input to \code{make_variants}.
//'
//' This function does NOT return an error if a sequence is chosen more times
//' than its length.
//'
//' @param total_mutations The total number of mutations (SNPs and indels).
//' @param seq_lens A vector of the sequence lengths.
//' @param seeds A vector seeds for the prng.
//'
//'
//' @return A numeric vector containing the number of mutations per sequence.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::vector<uint> sample_seqs(const uint& total_mutations,
                              const std::vector<double>& seq_lens,
                              const std::vector<uint>& seeds) {

    const uint n_seqs = seq_lens.size();
    const uint n_cores = seeds.size();

    // Creating output vector
    std::vector<uint> out_vec(n_seqs, 0);

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

    // Alias-sampling object
    const alias_FL FL(seq_lens);
    // sitmo prng
    sitmo::prng_engine eng(active_seed);

    // Parallelize the Loop
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint i = 0; i < total_mutations; i++){
        uint ind = alias_sample(n_seqs, FL, eng);
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





//' C++ equivalent of R's \code{choose} function.
//'
//' \emph{Note}: This function is not exported to R.
//'
//' @param n Integer value.
//' @param k Integer value.
//'
//' @return Binomial coefficient (also integer).
//'
//' @noRd
//'
inline uint cpp_choose(uint n, uint k) {
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




//' Calculate mean pairwise differences between samples using a vector of nucleotide
//'     frequencies.
//'
//'
//' @param sample_segr Vector of nucleotide frequencies at a given segregating site for
//'     all samples.
//'
//' @return Mean of the pairwise differences.
//'
//' @noRd
//'
double cpp_mean_pairwise_freqs(const std::vector<int>& sample_segr) {

    int N_nt = sample_segr.size();

    if (N_nt != 4) {
        stop("Vector input to cpp_mean_pairwise_freqs must be of length 4.");
    }

    int N = std::accumulate(sample_segr.begin(), sample_segr.end(), 0.0);

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
//' @noRd
//'
// [[Rcpp::export]]
List cpp_nt_freq(int N) {

    int a;
    std::vector<std::vector<int>> combos;
    std::vector<int> tmp_combos(4);
    std::vector<std::string> seq;
    std::vector<int> segr_bases(N);
    std::vector<double> mean_pws;
    double mean_pw;

    for (int d = 0; d <= floor(N / 4); d++) {
        for (int c = d; c <= floor((N - d) / 3); c++) {
            for (int b = c; b <= floor((N - c - d) / 2); b++) {
                a = N - b - c - d;
                if (a == N) continue;  // no reason to have all the same nucleotides
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


// These functions are used for iterating through sequences, ultimately creating a
// VarSet object.


// Change variant_set and return length for one segregating site

uint one_mutation(
        VarSet& var_set, const uint& seq_index,
        const sint& n, const uint& N, const uint& S, const sint64& current_pos,
        const uint& n_vars, const uint& seq_len,
        const std::vector<std::vector<uint>>& snp_combo_list,
        const std::vector<uint>& mutation_types,
        const std::vector<uint>& mutation_sizes,
        const alias_FL& fl, const uint& alias_n,
        sitmo::prng_engine& engine,
        const double& n2N = 50,
        const double& alpha = 0.8) {

    uint length;

    // Storing nucleos, sites, and sequence modifier
    std::string nucleos;
    std::vector<uint> sites;

    // Sampling for which type of mutation (SNP, insertion, deletion)
    uint mut_ind = alias_sample(alias_n, fl, engine);
    uint mut_type = mutation_types[mut_ind];

    // SNP
    if (mut_type == 0) {
        length = 0;
        std::vector<uint> combo = snp_combo_list[mut_ind];
        std::shuffle(combo.begin(), combo.end(), engine);
        // Creating and filling string
        nucleos = std::string(n_vars, 'x');
        for (uint j = 0, k = 0; j < 4; j++) {
            while (combo[j] > 0) {
                nucleos[k] = variants::bases[j];
                combo[j]--;
                k++;
            }
        }
        std::shuffle(nucleos.begin(), nucleos.end(), engine);
        for (uint v = 0; v < n_vars; v++) {
            VarSequence& vs(var_set[v][seq_index]);
            uint pos = static_cast<uint>(current_pos);
            if (pos >= vs.size()) continue;
            vs.add_substitution(nucleos[v], pos);
        }
    // InDel: Insertion or Deletion
    } else {
        // Make vector of the variants that have this indel:
        uint n_w_indel = ((double) engine() / variants::sitmo_max) * (n_vars - 1);
        std::vector<uint> w_indel(n_w_indel);
        vitter_d<std::vector<uint>>(w_indel, n_vars, engine, n2N, alpha);
        // Insertion
        if (mut_type == 1) {
            length = mutation_sizes[mut_ind];
            // Creating and filling string
            nucleos = std::string(length, 'x');
            for (uint j = 0; j < length; j++) {
                uint rnd = static_cast<double>(engine()) / variants::sitmo_max * 4;
                nucleos[j] = variants::bases[rnd];
            }
            for (uint v : w_indel) {
                VarSequence& vs(var_set[v][seq_index]);
                uint pos = static_cast<uint>(current_pos);
                if (pos >= vs.size()) continue;
                vs.add_insertion(nucleos, pos);
            }
            // Changing back to 0 bc insertions don't need to be skipped over like
            // deletions bc they don't take up existing nucleotides
            length = 0;
        // Deletion
        } else {
            length = mutation_sizes[mut_ind];
            // Make sure it doesn't span over the # positions left
            if (length > N - S - n + 1) length = N - S - n + 1;
            for (uint v : w_indel) {
                VarSequence& vs(var_set[v][seq_index]);
                uint pos = static_cast<uint>(current_pos);
                if (pos >= vs.size()) continue;
                vs.add_deletion(length, pos);
            }
            length--;
        }
    }
    return length;
}



//' Iterate and mutate one sequence.
//'
//' @noRd
//'
void one_seq(
        VarSet& var_set,
        const uint& seq_ind,
        const uint& n_muts,
        const uint& n_vars,
        const std::vector<std::vector<uint>>& snp_combo_list,
        const std::vector<uint> mutation_types,
        const std::vector<uint>& mutation_sizes,
        const alias_FL& fl,
        sitmo::prng_engine& engine,
        const double& n2N = 50,
        const double& alpha = 0.8
    ) {

    uint alias_n = fl.size();

    // These values are copied bc n and N will be changing
    sint n = n_muts;
    uint seq_len = var_set.reference[seq_ind].size();
    uint N = seq_len;

    // The # positions to skip before taking the next one (0 to (N - n - 1))
    uint S;
    // Keeping track of the current position
    sint64 current_pos = -1;  // (starts at -1 so it can reach 0 on the first skip)

    uint length; // Length of segregating sites

    // Stores function to sequentially find `S`
    std::function<uint(const sint&, const uint&, sitmo::prng_engine&,
                       const double)> algorithm;

    if (((n * n) / N) > n2N) {
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
        length = one_mutation(var_set, seq_ind, n, N, S, current_pos, n_vars, seq_len,
                              snp_combo_list, mutation_types, mutation_sizes,
                              fl, alias_n, engine, n2N, alpha);
        current_pos += length;
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
//'     for each sequence.
//' @param reference External pointer to a C++ \code{SequenceSet} object that
//'     represents the reference genome.
//' @param snp_combo_list Matrix of all possible nucleotide combinations among all
//'     variants per SNP.
//' @param snp_probs_cumsum Vector of sampling probabilities for each row in
//'     \code{snp_combo_list}.
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
//' @noRd
//'
//[[Rcpp::export]]
SEXP make_variants_(
        const std::vector<uint>& n_mutations,
        const SEXP& ref_xptr,
        const std::vector<std::vector<uint>>& snp_combo_list,
        const std::vector<double>& mutation_probs,
        const std::vector<uint>& mutation_types,
        const std::vector<uint>& mutation_sizes,
        std::vector<uint> seeds,
        double n2N = 50,
        double alpha = 0.8
    ) {

    const XPtr<RefGenome> reference(ref_xptr);

    const std::vector<uint> seq_lens = reference->seq_sizes();
    const uint n_seqs = reference->size();
    const uint n_cores = seeds.size();
    const uint n_vars = std::accumulate(snp_combo_list[0].begin(),
                                        snp_combo_list[0].end(), 0.0);

    if (n_mutations.size() != n_seqs) stop("n_mutations is incorrect length.");
    if (alpha > 1 || alpha < 0) stop("Invalid alpha. It must be [0,1].");

    XPtr<VarSet> var_set_xptr(new VarSet((*reference), n_vars), true);
    VarSet& var_set(*var_set_xptr);

    #ifdef _OPENMP
    #pragma omp parallel shared(var_set) num_threads(n_cores) if (n_cores > 1)
    {
    #endif

    uint active_seed;

    // Alias-sampling object
    const alias_FL fl(mutation_probs);

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
    for (uint s = 0; s < n_seqs; s++) {
        one_seq(var_set, s, n_mutations[s], n_vars,
                snp_combo_list, mutation_types, mutation_sizes, fl, engine,
                n2N, alpha);
    }
    #ifdef _OPENMP
    }
    #endif


    return var_set_xptr;
}


