#ifndef __JACKAL_TABLE_SAMPLER_H
#define __JACKAL_TABLE_SAMPLER_H


/*
 ********************************************************

 Table sampling from...
 Marsaglia, G., W. W. Tsang, and J. Wang. 2004. Fast generation of discrete random
 variables. Journal of Statistical Software 11.

 ********************************************************
 */

#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <pcg/pcg_random.hpp> // pcg prng

#include "jackal_types.h" // integer types
#include "pcg.h"  // pcg seeding
#include "util.h"  // decreasing_indices, str_stop


using namespace Rcpp;

namespace table_sampler {

    const std::string bases = "TCAG";

}



class TableSampler {
public:
    // Stores vectors of each category's Pr(sampled):
    std::vector<std::vector<uint32>> T = std::vector<std::vector<uint32>>(4);
    // Stores values at which to transition between vectors of `T`:
    std::vector<uint64> t = std::vector<uint64>(3, 0);
    /*
     For whether probs are evenly divisible into 2^16 bits, which means that only
     the first `T` vector should be used.
     This allows me to avoid using a 128-bit integer type.
    */
    bool even_probs = false;

    TableSampler() {};
    TableSampler(const std::vector<long double>& probs) {
        construct(probs);
    }
    TableSampler(const std::vector<double>& probs) {
        std::vector<long double> probs_;
        probs_.reserve(probs.size());
        for (uint32 i = 0; i < probs.size(); i++) {
            probs_.push_back(static_cast<long double>(probs[i]));
        }
        construct(probs_);
    }
    // Copy constructor
    TableSampler(const TableSampler& other) : T(other.T), t(other.t) {}

    uint32 sample(pcg64& eng) const;

private:

    uint32 dg(const uint64& m, const uint32& k) {
        uint64 x = ((m>>(64-16*k))&65535);
        return static_cast<uint32>(x);
    }

    /*
     Fill vector of integers that should sum to 2^64
     If only one outcome is possible (all but one p is zero), it makes `ints` have
     only one item, the index for that outcome
     */
    void fill_ints(const std::vector<long double>& p,
                          std::vector<uint64>& ints);

    // Most of the construction of the TableSampler object:
    void construct(const std::vector<long double>& probs);

};



/*
 EXAMPLE USAGE:

uint32 sample_rare_(SEXP xptr_sexp, const uint64& N, const uint32& rare) {

    XPtr<TableSampler> xptr(xptr_sexp);

    uint32 rares = 0;

    pcg32 eng = seeded_pcg();

    for (uint64 i = 0; i < N; i++) {
        uint32 k = xptr->sample(eng);
        if (k == rare) rares++;
    }

    return rares;
}

*/


/*
 Class template for table sampling a string, using an underlying TableSampler object.
 `chars_in` should be the characters to sample from, `probs` the probabilities of
 sampling those characters.
 `T` can be `std::string` or `RefSequence`. Others may work, but are not guaranteed.
 */
template <typename T>
class TableStringSampler {
public:

    T characters;

    TableStringSampler(const T& chars_in, const std::vector<double>& probs)
        : characters(chars_in), uint_sampler(probs), n(probs.size()) {
        if (probs.size() != chars_in.size()) {
            str_stop({"For a TableStringSampler construction, arguments probs and ",
                     "chars_in must be same length."});
        }
    }
    TableStringSampler() {}
    // copy constructor
    TableStringSampler(const TableStringSampler& other)
        : characters(other.characters), uint_sampler(other.uint_sampler),
          n(other.n) {}

    void sample(std::string& str, pcg64& eng) const {
        for (uint32 i = 0; i < str.size(); i++) {
            uint32 k = uint_sampler.sample(eng);
            str[i] = characters[k];
        }
        return;
    }
    char sample(pcg64& eng) const {
        uint32 k = uint_sampler.sample(eng);
        return characters[k];
    }

private:
    TableSampler uint_sampler;
    uint32 n;
};



#endif
