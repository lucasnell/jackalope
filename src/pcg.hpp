#ifndef __GEMINO_PCG_H
#define __GEMINO_PCG_H

#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include "pcg/pcg_extras.hpp"  // pcg 128-bit integer type
#include <pcg/pcg_random.hpp> // pcg prng

#include "gemino_types.hpp" // integer types


using namespace Rcpp;




namespace pcg {
    const double max32 = static_cast<double>(pcg32::max());
    const long double max64 = static_cast<long double>(pcg64::max());
}



/*
 ========================

 Seeding

 ========================
 */



// To sample for seeds before multi-core operations
inline std::vector<std::vector<uint64>> mc_seeds(const uint32& n_cores) {

    std::vector<std::vector<uint64>> sub_seeds(n_cores, std::vector<uint64>(8));

    for (uint32 i = 0; i < n_cores; i++) {
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
 -----------
 32-bit versions
 These require half as many 32-bit integers bc they're seeded from two 64-bit integers
 -----------
 */

// To sample for seeds before multi-core operations
inline std::vector<std::vector<uint64>> mc_seeds32(const uint32& n_cores) {

    std::vector<std::vector<uint64>> sub_seeds(n_cores, std::vector<uint64>(4));

    for (uint32 i = 0; i < n_cores; i++) {
        sub_seeds[i] = as<std::vector<uint64>>(Rcpp::runif(4,0,4294967296));
    }

    return sub_seeds;
}


// Fill two 64-bit seeds from four 32-bit seeds (casted to 64-bit)
inline void fill_seeds32(const std::vector<uint64>& sub_seeds,
                       uint64& seed1, uint64& seed2) {
    // Converting to two 64-bit seeds for pcg32
    seed1 = (sub_seeds[0]<<32) + sub_seeds[1];
    seed2 = (sub_seeds[2]<<32) + sub_seeds[3];
}

/*
 For single-core operations, you can use R's RNG for 32-bit random number generation.
 */
inline pcg32 seeded_pcg32() {

    // Four 32-bit seeds from unif_rand
    std::vector<uint64> sub_seeds = as<std::vector<uint64>>(Rcpp::runif(4,0,4294967296));

    uint64 seed1;
    uint64 seed2;
    fill_seeds32(sub_seeds, seed1, seed2);

    pcg32 out(seed1, seed2);
    return out;
}

/*
 For multi-core operations, you should call `mc_seeds` when outside multi-core mode,
 then input an inner `std::vector<uint32>` (from inside the object output from `mc_seeds`)
 to this function when in multi-core mode to seed the PRNG.
 */
// sub_seeds needs to be at least 4-long!
inline pcg32 seeded_pcg32(const std::vector<uint64>& sub_seeds) {

    uint64 seed1;
    uint64 seed2;
    fill_seeds32(sub_seeds, seed1, seed2);

    pcg32 out(seed1, seed2);
    return out;
}






/*
 ========================

 Number generation

 ========================
 */

// uniform in range [0,1]
inline double runif_0011(pcg32& eng) {
    return static_cast<double>(eng()) / pcg::max32;
}
inline long double runif_0011(pcg64& eng) {
    return static_cast<long double>(eng()) / pcg::max64;
}
// uniform in range [0,1)
inline double runif_001(pcg32& eng) {
    return static_cast<double>(eng()) / (pcg::max32 + 1);
}
inline long double runif_001(pcg64& eng) {
    return static_cast<long double>(eng()) / (pcg::max64 + 1);
}
// uniform in range (0,1)
inline double runif_01(pcg32& eng) {
    return (static_cast<double>(eng()) + 1) / (pcg::max32 + 2);
}
inline long double runif_01(pcg64& eng) {
    return (static_cast<long double>(eng()) + 1) / (pcg::max64 + 2);
}
// uniform in range (a,b)
inline double runif_ab(pcg32& eng, const double& a, const double& b) {
    return a + ((static_cast<double>(eng()) + 1) / (pcg::max32 + 2)) * (b - a);
}
inline long double runif_ab(pcg64& eng, const long double& a, const long double& b) {
    return a + ((static_cast<long double>(eng()) + 1) / (pcg::max64 + 2)) * (b - a);
}
// uniform in range [a,b]
inline uint32 runif_aabb(pcg32& eng, const uint32& a, const uint32& b) {
    return a + (static_cast<double>(eng()) / (pcg::max32 + 1)) * (b - a + 1);
}
inline uint64 runif_aabb(pcg64& eng, const uint64& a, const uint64& b) {
    return a + (static_cast<long double>(eng()) / (pcg::max64 + 1)) * (b - a + 1);
}
inline uint32 runif_aabb(pcg64& eng, const uint32& a, const uint32& b) {
    return a + (static_cast<long double>(eng()) / (pcg::max64 + 1)) * (b - a + 1);
}





#endif
