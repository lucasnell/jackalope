#ifndef __JACKAL_MUTATOR_H
#define __JACKAL_MUTATOR_H


/*
 Combining samplers for location and for mutation type into a mutation sampler for
 a single sequence.
 */


#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // gamma_distribution


#include "jackalope_types.h" // integer types
#include "seq_classes_var.h"  // Var* classes
#include "pcg.h"  // pcg seeding
#include "alias_sampler.h"  // alias method of sampling
#include "mutator_location.h"  // mutator location sampling classes
#include "mutator_type.h"      // mutator type sampling classes
#include "util.h"  // str_stop



using namespace Rcpp;





/*
 MutationSampler combines objects for sampling mutation types and new
 nucleotides for insertions.
 */
class MutationSampler {

    /*
     Create a new string of nucleotides (for insertions) of a given length and using
     an input rng engine
     */
    inline std::string new_nucleos(const uint64& len, pcg64& eng) const {
        std::string str(len, 'x');
        insert.sample(str, eng);
        return str;
    }
    // Does most of the work for the mutate methods (all but location sampling)
    inline double mutate__(pcg64& eng, const uint64& start, sint64& end);

public:

    // VarSequence object pointer to be manipulated
    VarSequence* var_seq;
    // For sampling the mutation location:
    LocationSampler location;
    // For sampling the type of mutation:
    MutationTypeSampler type;
    // For new insertion sequences:
    AliasStringSampler<std::string> insert;

    MutationSampler() {}

    MutationSampler(const MutationSampler& other)
        : var_seq(other.var_seq), location(other.location), type(other.type),
          insert(other.insert) {}

    MutationSampler& operator=(const MutationSampler& other) {
        if (other.var_seq) var_seq = other.var_seq;
        location = other.location;
        type = other.type;
        insert = other.insert;
        return *this;
    }

    void new_seq(VarSequence& vs_, const arma::mat& gamma_mat) {
        var_seq = &vs_;
        location.new_seq(vs_, gamma_mat);
        return;
    }

    // Add mutation and return the change in the sequence rate that results
    double mutate(pcg64& eng);

    /*
     Overloaded for only mutating within a range.
     It also updates `end` if an indel occurs in the range.
     Make sure to keep checking for situation where `end < start` (i.e., sequence section
     is empty).
     `// ***` mark difference between this and previous `mutate` versions
     */
    double mutate(pcg64& eng, const uint64& start, sint64& end);

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
