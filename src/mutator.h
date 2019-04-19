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
#include "site_var.h"  // SequenceGammas class
#include "weighted_reservoir.h"  // weighted_reservoir_* functions
#include "mutator_location.h"  // mutator location sampling classes
#include "mutator_type.h"      // mutator type sampling classes
#include "util.h"  // str_stop



using namespace Rcpp;





/*
 OneSeqMutationSampler combines objects for sampling mutation types and new
 nucleotides for insertions.

 Class `C` should be `LocationSampler` or `ChunkLocationSampler`.
 */
template <class C>
class OneSeqMutationSampler {

    /*
     Sample for mutation location based on rates by sequence region and nucleotide,
     for the whole sequence or a range.
     */
    inline uint32 sample_location(pcg64& eng,
                                  const uint32& start = 0, const uint32& end = 0,
                                  const bool& ranged = false) {
        return location.sample(eng, start, end, ranged);
    }

    /*
    Sample for mutation type based on nucleotide and rng engine
    */
    inline MutationInfo sample_type(const char& c, pcg64& eng) const {
        return type.sample(c, eng);
    }

    /*
    Create a new string of nucleotides (for insertions) of a given length and using
    an input rng engine
    */
    inline std::string new_nucleos(const uint32& len, pcg64& eng) const {
        std::string str(len, 'x');
        insert.sample(str, eng);
        return str;
    }

public:

    // VarSequence object pointer to be manipulated
    VarSequence* var_seq;
    // For sampling the mutation location:
    C location;
    // For sampling the type of mutation:
    MutationTypeSampler type;
    // For new insertion sequences:
    AliasStringSampler<std::string> insert;

    OneSeqMutationSampler() {}

    OneSeqMutationSampler(VarSequence& vs_,
                          const C& location_,
                          const MutationTypeSampler& type_,
                          const AliasStringSampler<std::string>& insert_)
        : var_seq(&vs_), location(location_), type(type_), insert(insert_) {}

    OneSeqMutationSampler(const OneSeqMutationSampler<C>& other)
        : var_seq(other.var_seq), location(other.location), type(other.type),
          insert(other.insert) {}

    OneSeqMutationSampler<C>& operator=(const OneSeqMutationSampler<C>& other) {
        if (other.var_seq) var_seq = other.var_seq;
        location = other.location;
        type = other.type;
        insert = other.insert;
        return *this;
    }

    void fill_ptrs(VarSequence& vs_) {
        var_seq = &vs_;
        location.fill_ptrs(vs_);
        return;
    }

    void fill_gamma(const arma::mat& gamma_mat) {
        location.mr().gammas = SequenceGammas(gamma_mat);
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
    double mutate(pcg64& eng, const uint32& start, sint64& end);



    double total_rate(const uint32& start = 0, const uint32& end = 0,
                      const bool& ranged = false) {
        return location.total_rate(start, end, ranged);
    }
};


// Shortening these names
typedef OneSeqMutationSampler<LocationSampler> MutationSampler;
typedef OneSeqMutationSampler<ChunkLocationSampler> ChunkMutationSampler;






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

    uint32 n_ins = insertion_rates.size();
    uint32 n_del = deletion_rates.size();
    uint32 n_muts = 4 + n_ins + n_del;

    // 1 vector of probabilities for each nucleotide: T, C, A, then G
    probs.resize(4);
    // Overall mutation rates for each nucleotide: T, C, A, then G
    q_tcag.reserve(4);

    for (uint32 i = 0; i < 4; i++) {

        std::vector<double>& qc(probs[i]);

        qc.reserve(n_muts);

        for (uint32 j = 0; j < Q.n_cols; j++) qc.push_back(Q(i, j));
        /*
         Make absolutely sure the diagonal is set to zero bc you don't want to
         mutate back to the same nucleotide
         */
        qc[i] = 0;

        // Add insertions, then deletions
        for (uint32 j = 0; j < n_ins; j++) {
            qc.push_back(insertion_rates[j] * 0.25);
        }
        for (uint32 j = 0; j < n_del; j++) {
            qc.push_back(deletion_rates[j] * 0.25);
        }
        // Get the overall mutation rate for this nucleotide
        double qi = std::accumulate(qc.begin(), qc.end(), 0.0);
        // Divide all in `qc` by `qi` to make them probabilities:
        for (uint32 j = 0; j < n_muts; j++) qc[j] /= qi;
        // Add `qi` to vector of rates by nucleotide:
        q_tcag.push_back(qi);
    }


    return;
}


//' Filling in mut_lengths vector
//'
//' @noRd
//'
inline void fill_mut_lengths(std::vector<sint32>& mut_lengths,
                             const std::vector<double>& insertion_rates,
                             const std::vector<double>& deletion_rates) {

    uint32 n_ins = insertion_rates.size();
    uint32 n_del = deletion_rates.size();
    uint32 n_muts = 4 + n_ins + n_del;

    // Now filling in mut_lengths vector
    mut_lengths.reserve(n_muts);
    for (uint32 i = 0; i < 4; i++) mut_lengths.push_back(0);
    for (uint32 i = 0; i < n_ins; i++) {
        mut_lengths.push_back(static_cast<sint32>(i+1));
    }
    for (uint32 i = 0; i < n_del; i++) {
        sint32 ds = static_cast<sint32>(i + 1);
        ds *= -1;
        mut_lengths.push_back(ds);
    }

    return;
}



//' Creates MutationSampler without any of the pointers.
//'
//'
//' `T` should be MutationSampler or ChunkMutationSampler
//' `T` should be LocationSampler or ChunkLocationSampler
//' MutationSampler should always go with LocationSampler, and
//' ChunkMutationSampler with ChunkLocationSampler
//'
//' Before actually using the object output from this function, make sure to...
//' * use `[Chunk]MutationSampler.fill_ptrs(VarSequence& var_seq)` to fill pointers.
//' * use `[Chunk]MutationSampler.fill_gamma(const arma::mat& gamma_mat)` to fill
//'   the gamma matrix.
//' * use `ChunkMutationSampler.location.change_chunk(chunk_size)` if using chunked
//'   version.
//'
//' @param Q A 4x4 matrix of substitution rates for each nucleotide.
//' @param pi_tcag Vector of nucleotide equilibrium frequencies for
//'     "T", "C", "A", and "G", respectively.
//' @param insertion_rates Vector of insertion rates.
//' @param deletion_rates Vector of deletion rates.
//'
//' @noRd
//'
template <typename T, typename U>
XPtr<T> make_mutation_sampler_base_(const arma::mat& Q,
                                    const std::vector<double>& pi_tcag,
                                    const std::vector<double>& insertion_rates,
                                    const std::vector<double>& deletion_rates) {

    std::vector<std::vector<double>> probs;
    std::vector<sint32> mut_lengths;
    std::vector<double> q_tcag;
    /*
    (1) Combine substitution, insertion, and deletion rates into a single vector
    (2) Fill the `q_tcag` vector with mutation rates for each nucleotide
    */
    fill_probs_q_tcag(probs, q_tcag, Q, pi_tcag, insertion_rates, deletion_rates);

    // Now filling in mut_lengths vector
    fill_mut_lengths(mut_lengths, insertion_rates, deletion_rates);

    /*
     Now create and fill output pointer to base sampler:
     */
    XPtr<T> out(new T());

    out->type = MutationTypeSampler(probs, mut_lengths);
    out->insert = AliasStringSampler<std::string>("TCAG", pi_tcag);

    MutationRates mr(q_tcag);
    out->location = U(mr);

    return out;
}








#endif
