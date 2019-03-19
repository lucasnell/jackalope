
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
#include "table_sampler.h"


using namespace Rcpp;



uint32 TableSampler::sample(pcg64& eng) const {
    uint64 j = eng();
    if (even_probs || j<t[0]) return T[0][j>>(64-16*1)];
    if (j<t[1]) return T[1][(j-t[0])>>(64-16*2)];
    if (j<t[2]) return T[2][(j-t[1])>>(64-16*3)];
    return T[3][j-t[2]];
}

void TableSampler::fill_ints(const std::vector<long double>& p,
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
void TableSampler::construct(const std::vector<long double>& probs) {

    uint32 n_tables = 4;

    uint32 n = probs.size();
    std::vector<uint64> ints;
    // Filling the `ints` vector based on `probs`
    this->fill_ints(probs, ints);

    // Taking care of scenario when just one output is possible
    if (ints.size() == 1) {

        // Below will result in sample only ever touching T[0]
        even_probs = true;
        /*
         Under this scenario, `fill_ints` makes `ints` only contain the index
         the we want returned.
         The lengths of T[0] and T[3] below will make sure that any value
         possible from a pcg64 RNG will return `ints[0]`
         */
        T[0] = std::vector<uint32>((1UL<<16), static_cast<uint32>(ints[0]));

    } else {

        std::vector<uint32> sizes(n_tables, 0);

        // Adding up sizes of `T` vectors:
        for (uint32 i = 0; i < n; i++) {
            for (uint32 k = 1; k <= n_tables; k++) {
                sizes[k-1] += this->dg(ints[i], k);
            }
        }

        // Dealing with situation where probabilities can be divided into 2^16 even bits:
        if (sizes[0] == 65536) {
            even_probs = true;
            // Remove unnecessary information:
            for (uint64& ii : ints) ii = dg(ii, 1);
            // Just fill first vector in `T`:
            T[0].reserve(65536);
            uint32 ind = 0; // index inside `T[k]`
            for (uint32 i = 0; i < n; i++) {
                const uint64& z(ints[i]);
                for (uint32 j = 0; j < z; j++) T[0].push_back(i);
                ind += z;
            }
        } else {
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
    }

    return;
}


