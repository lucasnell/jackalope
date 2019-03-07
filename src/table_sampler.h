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
#include "util.h"  // decreasing_indices


using namespace Rcpp;

namespace table_sampler {

    const std::string bases = "TCAG";


}



class TableSampler {
public:
    // Stores vectors of each category's Pr(sampled):
    std::vector<std::vector<uint32>> T;
    // Stores values at which to transition between vectors of `T`:
    std::vector<uint128> t;

    TableSampler() {};
    TableSampler(const std::vector<long double>& probs) : T(4), t(3, 0) {
        construct(probs);
    }
    TableSampler(const std::vector<double>& probs) : T(4), t(3, 0) {
        std::vector<long double> probs_;
        probs_.reserve(probs.size());
        for (uint32 i = 0; i < probs.size(); i++) {
            probs_.push_back(static_cast<long double>(probs[i]));
        }
        construct(probs_);
    }
    // Copy constructor
    TableSampler(const TableSampler& other) : T(other.T), t(other.t) {}

    inline uint32 sample(pcg64& eng) const {
        uint64 j = eng();
        if (j<t[0]) return T[0][j>>(64-16*1)];
        if (j<t[1]) return T[1][(j-t[0])>>(64-16*2)];
        if (j<t[2]) return T[2][(j-t[1])>>(64-16*3)];
        return T[3][j-t[2]];
    }

    inline uint32 sample(pcg32& eng) const {
        // Generate 64-bit integer from two draws of 32-bit RNG:
        uint64 j = eng();
        j <<= 32;
        j += eng();
        // Now sample as normal:
        if (j<t[0]) return T[0][j>>(64-16*1)];
        if (j<t[1]) return T[1][(j-t[0])>>(64-16*2)];
        if (j<t[2]) return T[2][(j-t[1])>>(64-16*3)];
        return T[3][j-t[2]];
    }

    void print() const {
        // names coincide with names from Marsaglia (2004)
        std::vector<std::string> names = {"AA", "BB", "CC", "DD"};
        for (uint32 i = 0; i < T.size(); i++) {
            Rcout << "T[" << i << "]:" << std::endl;
            for (const uint32& tt : T[i]) Rcout << tt << ' ';
            Rcout << std::endl;
        }
        Rcout << "t" << std::endl;
        for (const uint128& tt : t) Rcout << static_cast<double>(tt) << ' ';
        Rcout << std::endl;
    }

private:

    uint32 dg(const uint128& m, const uint32& k) {
        uint128 x = ((m>>(64ULL-16ULL*static_cast<uint128>(k)))&65535ULL);
        uint32 y = static_cast<uint32>(x);
        return y;
    }

    inline void fill_ints(const std::vector<long double>& p,
                          std::vector<uint128>& ints) {

        long double max_int = 18446744073709551616.0;  // 2^64

        uint32 n = p.size();

        std::vector<long double> pp(n);
        long double p_sum = std::accumulate(p.begin(), p.end(), 0.0);
        for (uint32 i = 0; i < n; i++) {
            long double x = max_int * p[i] / p_sum;
            pp[i] = std::round(x);
        }

        std::vector<uint32> inds = decreasing_indices<long double>(p);

        long double pp_sum = std::accumulate(pp.begin(), pp.end(), 0.0);
        long double d = max_int - pp_sum;

        uint32_t i = 0;
        // We need to remove from `pp`
        while (d < 0) {
            pp[inds[i]]--;
            i++;
            if (i == inds.size()) i = 0;
            d++;
        }
        // We need to add to `pp`
        while (d > 0) {
            pp[inds[i]]++;
            i++;
            if (i == inds.size()) i = 0;
            d--;
        }

        // Now fill `ints`:
        ints.reserve(n);
        for (uint32 i = 0; i < n; i++) {
            ints.push_back(static_cast<uint128>(pp[i]));
        }

        return;
    }

    // Most of the construction of the TableSampler object:
    inline void construct(const std::vector<long double>& probs) {

        uint32 n_tables = T.size();

        uint32 n = probs.size();
        std::vector<uint128> ints;
        // Filling the `ints` vector based on `probs`
        this->fill_ints(probs, ints);

        std::vector<uint32> sizes(n_tables, 0);
        // Adding up sizes of `T` vectors:
        for (uint32 i = 0; i < n; i++) {
            for (uint32 k = 1; k <= n_tables; k++) {
                sizes[k-1] += this->dg(ints[i], k);
            }
        }

        // Adding up thresholds in the `t` vector
        for (uint64 k = 0; k < (n_tables - 1); k++) {
            t[k] = sizes[k];
            t[k] <<= (64 - 16 * (1 + k));
            if (k > 0) t[k] += t[k-1];
        }

        // Taking care of scenario when just one output is possible
        if (std::accumulate(sizes.begin(), sizes.end(), 0ULL) == 0ULL) {
            // So it's always TRUE for the first `if ()` statement in sample:
            t[0] = (1ULL<<63);
            t[0] *= 2;
            /*
            Now filling in the index to the output with P = 1 so that sample
            always returns it. Bc we're iterating by 2^16, that's the number of
            items I have to fill in for the T[0] vector.
            */
            uint32 max_ind = 0;
            for (uint32 i = 0; i < ints.size(); i++) {
                if (ints[i] == t[0]) {
                    max_ind = i;
                    break;
                }
            }
            T[0] = std::vector<uint32>((1UL<<16), max_ind);
        } else {
            // Re-sizing `T` vectors:
            for (uint32 i = 0; i < n_tables; i++) T[i].resize(sizes[i]);
            // Filling `T` vectors
            for (uint32 k = 1; k <= n_tables; k++) {
                uint32 ind = 0; // index inside `T[k-1]`
                for (uint32 i = 0; i < n; i++) {
                    uint32 z = this->dg(ints[i], k);
                    for (uint32 j = 0; j < z; j++) T[k-1][ind + j] = i;
                    ind += z;
                }
            }
        }

        return;
    }
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
        for (uint32 i = 0; i < str.size(); i++) {
            uint32 k = uint_sampler.sample(eng);
            str[i] = characters[k];
        }
        return;
    }
    void sample(std::string& str, pcg64& eng) const {
        for (uint32 i = 0; i < str.size(); i++) {
            uint32 k = uint_sampler.sample(eng);
            str[i] = characters[k];
        }
        return;
    }

private:
    TableSampler uint_sampler;
    uint32 n;
};



#endif
