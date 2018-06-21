
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
#include "table_sampler.h"  // pcg seeding


using namespace Rcpp;



/*
 EXAMPLE USAGE: Sampling for a "rare" (low p) integer named `rare`

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
 Converts a `p` vector of probabilities to a `ints` vector of integers, where each
 integer represents the approximate expected value of "successes" from 2^32 runs.
 If the sum of `ints` is != 2^32, this randomly chooses values to change based on
 probabilities in `p`.
 I did it this way because larger probabilities should be less affected by changing
 their respective values in `ints`.
 */
inline void fill_ints(const std::vector<double>& p, std::vector<uint64>& ints,
                      pcg32& eng) {

    // Vector holding the transitory values that will eventually be inserted into `int`
    arma::vec pp(p);
    pp /= arma::accu(pp);
    pp *= static_cast<double>(1UL<<32);
    pp = arma::round(pp);

    // Converting to `ints`
    ints = arma::conv_to<std::vector<uint64>>::from(pp);

    double d = static_cast<double>(1UL<<32) - arma::accu(pp);

    // Vector for weighted sampling from vector of probabilities
    arma::vec p2(p);
    p2 /= arma::accu(p2);
    /*
     I'm not going to sample rare probabilities so anything < 2^-8 is set to zero
     for this sampling
     */
    double z = 1 / std::pow(2, 8);
    arma::uvec iv = arma::find(p2 < z);
    /*
     If there aren't any *above* this threshold, keep adding 8 to `x` in the expression
     `2^-x` until we would no longer be setting all probabilities to zero
     */
    while (iv.n_elem == p2.n_elem) {
        for (uint32 zz = 0; zz < 8; zz++) z /= 2;
        iv = arma::find(p2 < z);
    }
    p2(iv).fill(0);
    p2 /= arma::accu(p2);
    p2 = arma::cumsum(p2);

    // We need to remove from `ints`
    while (d < 0) {
        double u = static_cast<double>(eng()) / pcg::max;
        iv = arma::find(p2 >= u, 1);
        ints[iv(0)]--;
        d++;
    }
    // We need to add to `ints`
    while (d > 0) {
        double u = static_cast<double>(eng()) / pcg::max;
        iv = arma::find(p2 >= u, 1);
        ints[iv(0)]++;
        d--;
    }

    return;
}




TableSampler::TableSampler(const std::vector<double>& probs) : T(4), t(3, 0) {

    pcg32 eng = seeded_pcg();

    uint32 n_tables = T.size();

    uint32 n = probs.size();
    std::vector<uint64> ints(n);
    // Filling the `ints` vector based on `probs`
    fill_ints(probs, ints, eng);

    std::vector<uint32> sizes(n_tables, 0);
    // Adding up sizes of `T` vectors:
    for (uint32 i = 0; i < n; i++) {
        for (uint32 k = 1; k <= n_tables; k++) {
            sizes[k-1] += dg(ints[i], k);
        }
    }
    // Adding up thresholds in the `t` vector
    for (uint32 k = 0; k < (n_tables - 1); k++) {
        t[k] = sizes[k];
        t[k] <<= (32-8*(1+k));
        if (k > 0) t[k] += t[k-1];
    }
    // Taking care of scenario when just one output is possible
    if (std::accumulate(sizes.begin(), sizes.end(), 0ULL) == 0ULL) {
        // So it's always TRUE for the first `if ()` statement in sample:
        t[0] = (1ULL<<32);
        /*
         Now filling in the index to the output with P = 1 so that sample
         always returns it. Bc we're iterating by 2^8, that's the number of
         items I have to fill in for the T[0] vector.
         */
        uint32 max_ind = std::find(ints.begin(), ints.end(), t[0]) - ints.begin();
        T[0] = std::vector<uint32>((1UL<<8), max_ind);
    } else {
        // Re-sizing `T` vectors:
        for (uint32 i = 0; i < n_tables; i++) T[i].resize(sizes[i]);
        // Filling `T` vectors
        for (uint32 k = 1; k <= n_tables; k++) {
            uint32 ind = 0; // index inside `T[k-1]`
            for (uint32 i = 0; i < n; i++) {
                uint32 z = dg(ints[i], k);
                for (uint32 j = 0; j < z; j++) T[k-1][ind + j] = i;
                ind += z;
            }
        }
    }
}




void TableSampler::print() const {
    // names coincide with names from Marsaglia (2004)
    std::vector<std::string> names = {"AA", "BB", "CC", "DD"};
    for (uint32 i = 0; i < T.size(); i++) {
        Rcout << "T[" << i << "]:" << std::endl;
        for (const uint32& tt : T[i]) Rcout << tt << ' ';
        Rcout << std::endl;
    }
    Rcout << "t" << std::endl;
    for (const uint64& tt : t) Rcout << tt << ' ';
    Rcout << std::endl;
}

