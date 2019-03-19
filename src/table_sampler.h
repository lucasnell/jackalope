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
    std::vector<uint128> t = std::vector<uint128>(3, 0);

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
        uint128 j = static_cast<uint128>(eng());
        if (j<t[0]) {
            j >>= (64-16*1);
            uint32 jj = static_cast<uint32>(j);
            return T[0][jj];
        }
        if (j<t[1]) {
            j -= t[0];
            j >>= (64-16*2);
            uint32 jj = static_cast<uint32>(j);
            return T[1][jj];
        }
        if (j<t[2]) {
            j -= t[1];
            j >>= (64-16*3);
            uint32 jj = static_cast<uint32>(j);
            return T[2][jj];
        }
        j -= t[2];
        uint32 jj = static_cast<uint32>(j);
        return T[3][jj];
    }

private:

    uint32 dg(const uint128& m, const uint8& k) {
        // uint128 x = ((m>>(64-16*k))&65535);
        uint128 x = m;
        x >>= (64 - 16 * k);
        x &= static_cast<uint128>(65535);
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
            double ppi = static_cast<double>(pp[i]); // needed to pass appveyor tests
            ints.push_back(static_cast<uint128>(ppi));
        }

        return;
    }

    // Most of the construction of the TableSampler object:
    inline void construct(const std::vector<long double>& probs) {

        uint8 n_tables = 4;

        uint32 n = probs.size();
        std::vector<uint128> ints;
        // Filling the `ints` vector based on `probs`
        this->fill_ints(probs, ints);

        std::vector<uint32> sizes(n_tables, 0);
        // Adding up sizes of `T` vectors:
        for (uint32 i = 0; i < n; i++) {
            for (uint8 k = 1; k <= n_tables; k++) {
                sizes[k-1] += this->dg(ints[i], k);
            }
        }

        // Adding up thresholds in the `t` vector
        for (uint8 k = 0; k < (n_tables - 1); k++) {
            t[k] = sizes[k];
            t[k] <<= (64 - 16 * (1 + k));
            if (k > 0) t[k] += t[k-1];
        }

        // Taking care of scenario when just one output is possible
        if (std::accumulate(sizes.begin(), sizes.end(), 0ULL) == 0ULL) {
            // So it's always TRUE for the first `if ()` statement in sample:
            t[0] = static_cast<uint128>(18446744073709551616.0);
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
            for (uint8 i = 0; i < n_tables; i++) T[i].resize(sizes[i]);
            // Filling `T` vectors
            for (uint8 k = 1; k <= n_tables; k++) {
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
