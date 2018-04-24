#ifndef __GEMINO_TABLE_SAMPLER_H
#define __GEMINO_TABLE_SAMPLER_H


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

#include "gemino_types.h" // integer types
#include "pcg.h"  // pcg seeding


using namespace Rcpp;

namespace table_sampler {
    static const std::string bases = "ACGT";
}



class TableSampler {
public:
    // Stores vectors of each category's Pr(sampled):
    std::vector<std::vector<uint>> T;
    // Stores values at which to transitiion between vectors of `T`:
    std::vector<uint> t;

    TableSampler();
    TableSampler(const std::vector<double>& probs);
    // Copy constructor
    TableSampler(const TableSampler& other) : T(other.T), t(other.t) {}

    uint sample(pcg32& eng) const;

    void print() const;

private:
    static uint dg(const uint& m, const uint& k) {
        uint x = ((m>>(32-8*k))&255);
        return x;
    }
};



/*
 EXAMPLE USAGE:

uint sample_rare_(SEXP xptr_sexp, const uint64& N, const uint& rare) {

    XPtr<TableSampler> xptr(xptr_sexp);

    uint rares = 0;

    pcg32 eng = seeded_pcg();

    for (uint64 i = 0; i < N; i++) {
        uint k = xptr->sample(eng);
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
            stop("For a TableStringSampler construction, arguments probs and chars_in ",
                 "must be same length.");
        }
    }
    TableStringSampler() {}
    // copy constructor
    TableStringSampler(const TableStringSampler& other)
        : characters(other.characters), uint_sampler(other.uint_sampler),
          n(other.n) {}

    void sample(std::string& str, pcg32& eng) const {
        for (uint i = 0; i < str.size(); i++) {
            uint k = uint_sampler.sample(eng);
            str[i] = characters[k];
        }
        return;
    }

private:
    TableSampler uint_sampler;
    uint n;
};



#endif
