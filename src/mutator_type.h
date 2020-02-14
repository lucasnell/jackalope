#ifndef __JACKALOPE_MUTATOR_TYPE_H
#define __JACKALOPE_MUTATOR_TYPE_H



/*
 This defines classes for sampling mutation type based on the starting
 nucleotide.

 Used only for segregating sites method.
 */

#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // gamma_distribution


#include "jackalope_types.h" // integer types
#include "hap_classes.h"  // Hap* classes
#include "pcg.h"  // pcg seeding
#include "alias_sampler.h"  // alias method of sampling
#include "util.h"  // str_stop



using namespace Rcpp;




/*
 For a single mutation's info
 */
struct MutationInfo {
    char nucleo;
    sint64 length;

    MutationInfo() : nucleo('X'), length() {}
    MutationInfo(const MutationInfo& other)
        : nucleo(other.nucleo), length(other.length) {}
    MutationInfo& operator=(const MutationInfo& other) {
        nucleo = other.nucleo;
        length = other.length;
        return *this;
    }

    // Initialize from an index and mut-lengths vector
    MutationInfo (const uint64& ind, const std::vector<sint64>& mut_lengths)
        : nucleo('\0'), length(0) {
        if (ind < 4) {
            nucleo = jlp::bases[ind];
        } else {
            length = mut_lengths[ind];
        }
    }
};



/*
 For constructors, this creates a vector of indices for each char in "TCAG" (0 to 3).
 Because char objects can be easily cast to uints, I can input a char from a chromosome
 and get out an index to which AliasSampler object to sample from.
 This way is much faster than using an unordered_map.
 Using 8-bit uints bc the char should never be >= 256.
 */
inline std::vector<uint8> make_base_inds() {
    std::vector<uint8> base_inds(256, 4);
    uint8 i = 0;
    for (const char& c : jlp::bases) {
        base_inds[c] = i;
        i++;
    }
    return base_inds;
}



/*
 For alias-sampling mutation types depending on which nucleotide you start with.
 The `mut_lengths` vector tells how long each mutation is.
 This field is 0 for substitions, < 0 for deletions, and > 0 for insertions.
 The `base_inds` field allows me to convert the characters 'T', 'C', 'A', or 'G'
 (cast to uints) into uints from 0 to 3.
 */
class MutationTypeSampler {

    std::vector<AliasSampler> sampler;
    std::vector<sint64> mut_lengths;
    std::vector<uint8> base_inds;

public:

    MutationTypeSampler() : sampler(4), mut_lengths(), base_inds(make_base_inds()) {};
    MutationTypeSampler(const std::vector<std::vector<double>>& probs,
                        const std::vector<sint64>& mut_lengths_)
    : sampler(4), mut_lengths(mut_lengths_), base_inds(make_base_inds()) {
        if (probs.size() != 4) stop("probs must be size 4.");
        for (uint64 i = 0; i < 4; i++) sampler[i] = AliasSampler(probs[i]);
    }
    // copy constructor
    MutationTypeSampler(const MutationTypeSampler& other)
        : sampler(other.sampler), mut_lengths(other.mut_lengths),
          base_inds(other.base_inds) {}
    // Assignment operator
    MutationTypeSampler& operator=(const MutationTypeSampler& other) {
        sampler = other.sampler;
        mut_lengths = other.mut_lengths;
        base_inds = make_base_inds();
        return *this;
    }

    /*
     Sample a mutation based on an input nucleotide.
     `c` gets cast to an uint64, which is then input to `base_inds` to get the index
     from 0 to 3.
     */
    MutationInfo sample(const char& c, pcg64& eng) const {
        uint8 j = base_inds[c];
        if (j > 3) return MutationInfo(); // only mutating T, C, A, or G
        uint64 ind = sampler[j].sample(eng);
        MutationInfo mi(ind, mut_lengths);
        return mi;
    }

};





//' Fill probs and q_tcag vectors.
//'
//' (1) Combine substitution, insertion, and deletion rates into a single vector
//' (2) Fill the `q_tcag` vector with mutation rates for each nucleotide
//'
//' @noRd
//'
inline void fill_probs_q_tcag(std::vector<std::vector<double>>& probs,
                              std::vector<double>& q_tcag,
                              const arma::mat& Q,
                              const std::vector<double>& pi_tcag,
                              const std::vector<double>& insertion_rates,
                              const std::vector<double>& deletion_rates) {

    uint64 n_ins = insertion_rates.size();
    uint64 n_del = deletion_rates.size();
    uint64 n_muts = 4 + n_ins + n_del;

    // 1 vector of probabilities for each nucleotide: T, C, A, then G
    probs.resize(4);
    // Overall mutation rates for each nucleotide: T, C, A, then G
    q_tcag.reserve(4);

    for (uint64 i = 0; i < 4; i++) {

        std::vector<double>& qc(probs[i]);

        qc.reserve(n_muts);

        for (uint64 j = 0; j < Q.n_cols; j++) qc.push_back(Q(i, j));
        /*
         Make absolutely sure the diagonal is set to zero bc you don't want to
         mutate back to the same nucleotide
         */
        qc[i] = 0;

        // Add insertions, then deletions
        for (uint64 j = 0; j < n_ins; j++) {
            qc.push_back(insertion_rates[j] * 0.25);
        }
        for (uint64 j = 0; j < n_del; j++) {
            qc.push_back(deletion_rates[j] * 0.25);
        }
        // Get the overall mutation rate for this nucleotide
        double qi = std::accumulate(qc.begin(), qc.end(), 0.0);
        // Divide all in `qc` by `qi` to make them probabilities:
        for (uint64 j = 0; j < n_muts; j++) qc[j] /= qi;
        // Add `qi` to vector of rates by nucleotide:
        q_tcag.push_back(qi);
    }


    return;
}


//' Filling in mut_lengths vector
//'
//' @noRd
//'
inline void fill_mut_lengths(std::vector<sint64>& mut_lengths,
                             const std::vector<double>& insertion_rates,
                             const std::vector<double>& deletion_rates) {

    uint64 n_ins = insertion_rates.size();
    uint64 n_del = deletion_rates.size();
    uint64 n_muts = 4 + n_ins + n_del;

    // Now filling in mut_lengths vector
    mut_lengths.reserve(n_muts);
    for (uint64 i = 0; i < 4; i++) mut_lengths.push_back(0);
    for (uint64 i = 0; i < n_ins; i++) {
        mut_lengths.push_back(static_cast<sint64>(i+1));
    }
    for (uint64 i = 0; i < n_del; i++) {
        sint64 ds = static_cast<sint64>(i + 1);
        ds *= -1;
        mut_lengths.push_back(ds);
    }

    return;
}









#endif
