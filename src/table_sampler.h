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





class TableSampler {
public:
    // Stores vectors of each category's Pr(sampled):
    std::vector<std::vector<uint>> T;
    // Stores values at which to transitiion between vectors of `T`:
    std::vector<uint> t;

    TableSampler(const std::vector<double>& probs, pcg32& eng);

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



#endif
