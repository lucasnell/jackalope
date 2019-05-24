#ifndef __JACKAL_MUTATOR_TYPE_H
#define __JACKAL_MUTATOR_TYPE_H



/*
 This defines classes for sampling mutation type based on the starting
 nucleotide.
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
#include "util.h"  // str_stop



using namespace Rcpp;

namespace mut_type {
    const std::string bases = "TCAG";
}




/*
 For a single mutation's info
 */
struct MutationInfo {
    char nucleo;
    sint64 length;

    MutationInfo() : nucleo(), length() {}
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
            nucleo = mut_type::bases[ind];
        } else {
            length = mut_lengths[ind];
        }
    }
};



/*
 For constructors, this creates a vector of indices for each char in "TCAG" (0 to 3).
 Because char objects can be easily cast to uints, I can input a char from a sequence
 and get out an index to which AliasSampler object to sample from.
 This way is much faster than using an unordered_map.
 Using 8-bit uints bc the char should never be >= 256.
 It's only of length 85 (versus 256 in LocationSampler class) because characters other
 than T, C, A, or G have rates hard-coded to zero and should never be chosen.
 */
inline std::vector<uint8> make_base_inds() {
    std::vector<uint8> base_inds(85);
    uint8 i = 0;
    for (const char& c : mut_type::bases) {
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
        uint64 ind = sampler[base_inds[c]].sample(eng);
        MutationInfo mi(ind, mut_lengths);
        return mi;
    }

};







#endif
