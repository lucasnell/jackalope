#ifndef __GEMINO_PCG_H
#define __GEMINO_PCG_H

#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <pcg/pcg_random.hpp> // pcg prng

#include "gemino_types.h" // integer types
#include "pcg.h" // pcg seeder


using namespace Rcpp;


namespace pcg {
    static const double max = static_cast<double>(pcg32::max());
}




// To sample for seeds before multi-core operations
inline std::vector<std::vector<uint>> mc_seeds(const uint& n_cores) {

    std::vector<std::vector<uint>> sub_seeds(n_cores, std::vector<uint>(4));

    for (uint i = 0; i < n_cores; i++) {
        sub_seeds[i] = as<std::vector<uint>>(Rcpp::runif(4,0,4294967296));
    }

    return sub_seeds;
}


/*
 For single-core operations, you can use R's RNG for 32-bit random number generation.
 */
inline pcg32 seeded_pcg() {

    // Four 32-bit seeds from unif_rand
    std::vector<uint64> sub_seeds = as<std::vector<uint64>>(Rcpp::runif(4,0,4294967296));

    // Converting to two 64-bit seeds for pcg32
    uint64 seed1 = (sub_seeds[0]<<32) + sub_seeds[1];
    uint64 seed2 = (sub_seeds[2]<<32) + sub_seeds[3];

    pcg32 out(seed1, seed2);
    return out;
}

/*
 For multi-core operations, you should call `mc_seeds` when outside multi-core mode,
 then input an inner `std::vector<uint>` (from inside the object output from `mc_seeds`)
 to this function when in multi-core mode to seed the PRNG.
 */
inline pcg32 seeded_pcg(const std::vector<uint>& sub_seeds) {

    // 32-bit seeds from input seed vector
    uint64 seed1 = (static_cast<uint64>(sub_seeds[0])<<32) + sub_seeds[1];
    uint64 seed2 = (static_cast<uint64>(sub_seeds[2])<<32) + sub_seeds[3];

    pcg32 out(seed1, seed2);
    return out;
}


#endif
