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
#include "util.h"  // decreasing_indices


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

    inline uint32 sample(pcg64& eng) const {
        uint64 j = eng();
        if (j<t[0]) return T[0][j>>(64-16*1)];
        if (j<t[1]) return T[1][(j-t[0])>>(64-16*2)];
        if (j<t[2]) return T[2][(j-t[1])>>(64-16*3)];
        return T[3][j-t[2]];
    }

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
    inline void fill_ints(const std::vector<long double>& p,
                          std::vector<uint64>& ints) {

        long double max_int = 18446744073709551616.0;  // 2^64

        uint32 n = p.size();

        std::vector<long double> pp(n);
        long double p_sum = std::accumulate(p.begin(), p.end(), 0.0);
        for (uint32 i = 0; i < n; i++) {
            long double x = p[i] / p_sum;
            if (x == 1) {
                ints = { static_cast<uint64>(i) };
                return;
            }
            x *= max_int;
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
            double ppi = static_cast<double>(pp[i]); // needed to pass appveyor tests
            ints.push_back(static_cast<uint64>(ppi));
        }

        return;
    }

    // Most of the construction of the TableSampler object:
    inline void construct(const std::vector<long double>& probs) {

        uint32 n_tables = 4;

        uint32 n = probs.size();
        std::vector<uint64> ints;
        // Filling the `ints` vector based on `probs`
        this->fill_ints(probs, ints);

        // Taking care of scenario when just one output is possible
        if (ints.size() == 1) {

            // Below will result in sample only ever touching T[0] and T[3]
            t = std::vector<uint64>(4, 18446744073709551615ULL);
            /*
             Under this scenario, `fill_ints` makes `ints` only contain the index
             the we want returned.
             The lengths of T[0] and T[3] below will make sure that any value
             possible from a pcg64 RNG will return `ints[0]`
             */
            T[0] = std::vector<uint32>((1UL<<16), static_cast<uint32>(ints[0]));
            T[3] = std::vector<uint32>(1, static_cast<uint32>(ints[0]));

        } else {

            std::vector<uint32> sizes(n_tables, 0);
            // Adding up sizes of `T` vectors:
            for (uint32 i = 0; i < n; i++) {
                for (uint32 k = 1; k <= n_tables; k++) {
                    sizes[k-1] += this->dg(ints[i], k);
                }
            }
            // Adding up thresholds in the `t` vector
            for (uint32 k = 0; k < (n_tables - 1); k++) {
                t[k] = sizes[k];
                t[k] <<= (64 - 16 * (1 + k));
                if (k > 0) t[k] += t[k-1];
            }
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
