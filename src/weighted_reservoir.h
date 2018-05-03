#ifndef __GEMINO_WEIGHTED_RESERVOIR_H
#define __GEMINO_WEIGHTED_RESERVOIR_H

/*
 ************************************************************************
 ************************************************************************

 Header containing templates to use for weighted reservoir sampling.

 Method is from...
    Efraimidis, P. S., and P. G. Spirakis. 2006. Weighted random sampling with a
    reservoir. Information Processing Letters 97:181–185.

 ************************************************************************
 ************************************************************************
 */


#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng

#include "gemino_types.h" // integer types
#include "pcg.h" // runif_* functions


/*
 ================================================================
 ================================================================

 Weighted reservoir sampling

 ================================================================
 ================================================================
 */

/*
 Templates that do weighted reservoir sampling for a section of a vector, returning
 one unsigned integer index to the sampled location in the `rates` vector.

 The first version keeps track of the number of basepairs that are being surveyed
 as it passes along the `rates` vector.
 It stops when the number of basepairs passed is greater than the `chunk_size` argument.
 This version also starts from the beginning of `rates` if it reaches the end but
 hasn't reached the `chunk_size` threshold.
 The number of basepairs is determined using the `size` method in class `T`,
 which takes an integer index as an argument and returns the size (in basepairs) of the
 region the index points to.

 Both versions require that class `T` has a bracket operator (`operator[]`) defined
 that returns the rate for a given index.
 They also both require that `end` is <= the maximum index for `rates`.

 */


template <typename T>
inline uint weighted_reservoir_chunk_(const uint& start, const uint& end,
                                      const uint& chunk_size,
                                      const T& rates,
                                      pcg32& eng) {

    double r, key, X, w, t;

    // Create objects to store currently-selected position and key
    r = runif_01(eng);
    key = std::log(r) / rates[start]; // log(key)
    double largest_key = key;  // largest key (the one we're going to keep)
    uint largest_pos = start;  // position where largest key was found

    uint c = start;
    uint n_bp = rates.size(c);
    while (n_bp < chunk_size) {
        r = runif_01(eng);
        X = std::log(r) / largest_key;  // largest_key is already logged
        uint i = c + 1;
        if (i > end) i = 0;
        double wt_sum0 = rates[c];
        double wt_sum1 = wt_sum0 + rates[i];
        n_bp += rates.size(i);
        while (X > wt_sum1 && n_bp < chunk_size) {
            i++;
            wt_sum0 += rates[(i-1)];
            if (i > end) i = 0;
            n_bp += rates.size(i);
            wt_sum1 += rates[i];
        }
        if (X > wt_sum1) break;
        if (wt_sum0 >= X) continue;

        largest_pos = i;

        w = rates[i];
        t = std::exp(w * largest_key); // key is log(key)
        r = runif_ab(eng, t, 1.0);
        key = std::log(r) / w; // log(key)
        largest_key = key;

        c = i;
    }

    return largest_pos;
}


template <typename T>
inline uint weighted_reservoir_(const uint& start, const uint& end,
                                const T& rates, pcg32& eng) {

    double r, key, X, w, t;


    // Create objects to store currently-selected position and key
    r = runif_01(eng);
    key = std::log(r) / rates[start]; // log(key)  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    double largest_key = key;  // largest key (the one we're going to keep)
    uint largest_pos = start;  // position where largest key was found

    uint c = start;
    while (c < end) {
        r = runif_01(eng);
        X = std::log(r) / largest_key;  // log(key)  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        uint i = c + 1;
        double wt_sum0 = rates[c];
        double wt_sum1 = wt_sum0 + rates[i];
        while (X > wt_sum1 && i < end) {
            i++;
            wt_sum0 += rates[(i-1)];
            wt_sum1 += rates[i];
        }
        if (X > wt_sum1) break;
        if (wt_sum0 >= X) continue;

        largest_pos = i;

        w = rates[i];
        t = std::exp(w * largest_key); // key is log(key)  // <<<<<<<<<<<<<<<<<<<<<<<<<<<
        r = runif_ab(eng, t, 1.0);
        key = std::log(r) / w; // log(key)  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        largest_key = key;

        c = i;
    }

    return largest_pos;
}


#endif
