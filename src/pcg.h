#ifndef __JACKALOPE_PCG_H
#define __JACKALOPE_PCG_H

#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include "pcg/pcg_extras.hpp"  // pcg 128-bit integer type
#include <pcg/pcg_random.hpp> // pcg prng

#include "jackalope_types.h" // integer types


using namespace Rcpp;




namespace pcg {
    const long double max64 = static_cast<long double>(pcg64::max());
}



/*
 ========================

 Seeding

 ========================
 */



// To sample for seeds before multi-thread operations
inline std::vector<std::vector<uint64>> mt_seeds(const uint64& n_threads) {

    std::vector<std::vector<uint64>> sub_seeds(n_threads, std::vector<uint64>(8));

    for (uint64 i = 0; i < n_threads; i++) {
        sub_seeds[i] = as<std::vector<uint64>>(Rcpp::runif(8,0,4294967296));
    }

    return sub_seeds;
}

inline void fill_seeds(const std::vector<uint64>& sub_seeds,
                         uint128& seed1, uint128& seed2) {

    uint128 seed64_1 = (sub_seeds[0]<<32) + sub_seeds[1];
    uint128 seed64_2 = (sub_seeds[2]<<32) + sub_seeds[3];
    uint128 seed64_3 = (sub_seeds[4]<<32) + sub_seeds[5];
    uint128 seed64_4 = (sub_seeds[6]<<32) + sub_seeds[7];

    seed1 = (seed64_1<<64) + seed64_2;
    seed2 = (seed64_3<<64) + seed64_4;

    return;
}

inline pcg64 seeded_pcg() {

    // 32-bit seeds from unif_rand
    std::vector<uint64> sub_seeds = as<std::vector<uint64>>(Rcpp::runif(8,0,4294967296));
    uint128 seed1;
    uint128 seed2;

    fill_seeds(sub_seeds, seed1, seed2);

    pcg64 out(seed1, seed2);
    return out;
}

// sub_seeds needs to be at least 8-long!
inline pcg64 seeded_pcg(const std::vector<uint64>& sub_seeds) {

    uint128 seed1;
    uint128 seed2;
    fill_seeds(sub_seeds, seed1, seed2);

    pcg64 out(seed1, seed2);

    return out;
}





/*
 ========================

 Number generation

 ========================
 */
// uniform in range (0,1)
inline long double runif_01(pcg64& eng) {
    return (static_cast<long double>(eng()) + 1) / (pcg::max64 + 2);
}
// uniform in range (a,b)
inline long double runif_ab(pcg64& eng, const long double& a, const long double& b) {
    return a + ((static_cast<long double>(eng()) + 1) / (pcg::max64 + 2)) * (b - a);
}





#endif
