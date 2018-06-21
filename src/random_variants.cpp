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
#include <pcg/pcg_random.hpp> // pcg prng
#ifdef _OPENMP
#include <omp.h>  // omp
#endif



#include "gemino_types.h"  // integer types
#include "sequence_classes.h"  // Ref* and Var* classes
#include "vitter_algorithms.h"  // vitter_d
#include "table_sampler.h" // table sampling
#include "pcg.h" // pcg sampler types
#include "util.h" // cpp_choose


using namespace Rcpp;







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
std::vector<uint32> sample_seqs(const uint32& total_mutations,
                              const std::vector<double>& seq_lens,
                              const uint32& n_cores) {

    const uint32 n_seqs = seq_lens.size();

    // Creating output vector
    std::vector<uint32> out_vec(n_seqs, 0);

    // Seeds
    const std::vector<std::vector<uint64>> seeds = mc_seeds(n_cores);

    // Table-sampling object
    const TableSampler sampler(seq_lens);

    #ifdef _OPENMP
    #pragma omp parallel default(shared) num_threads(n_cores) if(n_cores > 1)
    {
    #endif

    std::vector<uint64> active_seeds;

    // Write the active seed per core or just write one of the seeds.
    #ifdef _OPENMP
    uint32 active_thread = omp_get_thread_num();
    active_seeds = seeds[active_thread];
    #else
    active_seeds = seeds[0];
    #endif

    // pcg prng
    pcg32 eng = seeded_pcg(active_seeds);

    // Parallelize the Loop
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint32 i = 0; i < total_mutations; i++){
        uint32 ind = sampler.sample(eng);
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

    uint32 N_nt = sample_segr.size();

    if (N_nt != 4) {
        stop("Vector input to cpp_mean_pairwise_freqs must be of length 4.");
    }

    int N = std::accumulate(sample_segr.begin(), sample_segr.end(), 0);

    double total_pairs = cpp_choose(static_cast<double>(N), 2.0);
    double same_pairs = 0;
    for (uint32 i = 0; i < N_nt; i++) {
        same_pairs += cpp_choose(static_cast<double>(sample_segr[i]), 2.0);
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

uint32 one_mutation(
        VarSet& var_set, const uint32& seq_index,
        const sint32& n, const uint32& N, const uint32& S, const sint64& current_pos,
        const uint32& n_vars, const uint32& seq_len,
        const std::vector<std::vector<uint32>>& snp_combo_list,
        const std::vector<uint32>& mutation_types,
        const std::vector<uint32>& mutation_sizes,
        const TableSampler& sampler,
        pcg32& engine,
        const double& n2N = 50,
        const double& alpha = 0.8) {

    uint32 length;

    // Storing nucleos, sites, and sequence modifier
    std::string nucleos;
    std::vector<uint32> sites;

    // Sampling for which type of mutation (SNP, insertion, deletion)
    uint32 mut_ind = sampler.sample(engine);
    uint32 mut_type = mutation_types[mut_ind];

    // SNP
    if (mut_type == 0) {
        length = 0;
        std::vector<uint32> combo = snp_combo_list[mut_ind];
        std::shuffle(combo.begin(), combo.end(), engine);
        // Creating and filling string
        nucleos = std::string(n_vars, 'x');
        for (uint32 j = 0, k = 0; j < 4; j++) {
            while (combo[j] > 0) {
                nucleos[k] = table_sampler::bases[j];
                combo[j]--;
                k++;
            }
        }
        std::shuffle(nucleos.begin(), nucleos.end(), engine);
        for (uint32 v = 0; v < n_vars; v++) {
            VarSequence& vs(var_set[v][seq_index]);
            uint32 pos = static_cast<uint32>(current_pos);
            if (pos >= vs.size()) continue;
            vs.add_substitution(nucleos[v], pos);
        }
    // InDel: Insertion or Deletion
    } else {
        // Make vector of the variants that have this indel:
        uint32 n_w_indel = runif_01(engine) * (n_vars - 1);
        std::vector<uint32> w_indel(n_w_indel);
        vitter_d<std::vector<uint32>>(w_indel, n_vars, engine, n2N, alpha);
        // Insertion
        if (mut_type == 1) {
            length = mutation_sizes[mut_ind];
            // Creating and filling string
            nucleos = std::string(length, 'x');
            for (uint32 j = 0; j < length; j++) {
                uint32 rnd = runif_01(engine) * 4;
                nucleos[j] = table_sampler::bases[rnd];
            }
            for (uint32 v : w_indel) {
                VarSequence& vs(var_set[v][seq_index]);
                uint32 pos = static_cast<uint32>(current_pos);
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
            for (uint32 v : w_indel) {
                VarSequence& vs(var_set[v][seq_index]);
                uint32 pos = static_cast<uint32>(current_pos);
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
        const uint32& seq_ind,
        const uint32& n_muts,
        const uint32& n_vars,
        const std::vector<std::vector<uint32>>& snp_combo_list,
        const std::vector<uint32> mutation_types,
        const std::vector<uint32>& mutation_sizes,
        const TableSampler& sampler,
        pcg32& engine,
        const double& n2N = 50,
        const double& alpha = 0.8
    ) {

    // These values are copied bc n and N will be changing
    sint32 n = n_muts;
    uint32 seq_len = var_set.reference[seq_ind].size();
    uint32 N = seq_len;

    // The # positions to skip before taking the next one (0 to (N - n - 1))
    uint32 S;
    // Keeping track of the current position
    sint64 current_pos = -1;  // (starts at -1 so it can reach 0 on the first skip)

    uint32 length; // Length of segregating sites

    // Stores function to sequentially find `S`
    std::function<uint32(const sint32&, const uint32&, pcg32&,
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
            S = runif_01(engine) * N;
        }
        current_pos += S + 1;
        // This function returns the length, but modifies variant_set
        length = one_mutation(var_set, seq_ind, n, N, S, current_pos, n_vars, seq_len,
                              snp_combo_list, mutation_types, mutation_sizes,
                              sampler, engine, n2N, alpha);
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
        const std::vector<uint32>& n_mutations,
        const SEXP& ref_xptr,
        const std::vector<std::vector<uint32>>& snp_combo_list,
        const std::vector<double>& mutation_probs,
        const std::vector<uint32>& mutation_types,
        const std::vector<uint32>& mutation_sizes,
        const uint32& n_cores,
        double n2N = 50,
        double alpha = 0.8
    ) {

    const XPtr<RefGenome> reference(ref_xptr);

    const std::vector<uint32> seq_lens = reference->seq_sizes();
    const uint32 n_seqs = reference->size();
    const std::vector<std::vector<uint64>> seeds = mc_seeds(n_cores);
    const uint32 n_vars = std::accumulate(snp_combo_list[0].begin(),
                                        snp_combo_list[0].end(), 0.0);

    if (n_mutations.size() != n_seqs) stop("n_mutations is incorrect length.");
    if (alpha > 1 || alpha < 0) stop("Invalid alpha. It must be [0,1].");

    XPtr<VarSet> var_set_xptr(new VarSet((*reference), n_vars), true);
    VarSet& var_set(*var_set_xptr);

    // Table-sampling object
    const TableSampler sampler(mutation_probs);

    #ifdef _OPENMP
    #pragma omp parallel shared(var_set) num_threads(n_cores) if (n_cores > 1)
    {
    #endif

    std::vector<uint64> active_seeds;

    // Write the active seed per core or just write one of the seeds.
    #ifdef _OPENMP
    active_seeds = seeds[omp_get_thread_num()];
    #else
    active_seeds = seeds[0];
    #endif

    pcg32 engine = seeded_pcg(active_seeds);

    // Parallelize the Loop
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint32 s = 0; s < n_seqs; s++) {
        one_seq(var_set, s, n_mutations[s], n_vars,
                snp_combo_list, mutation_types, mutation_sizes, sampler, engine,
                n2N, alpha);
    }
    #ifdef _OPENMP
    }
    #endif


    return var_set_xptr;
}


