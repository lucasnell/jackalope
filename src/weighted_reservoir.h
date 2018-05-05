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
#include "vitter_algorithms.h" // vitter_d


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

 NOTE: Make sure `distr` is constructed with rate = 1, like the following:
    `std::exponential_distribution<double> distr(1.0);`
    It should never have a different rate!
 */




template <typename T>
inline uint weighted_reservoir_chunk_(const uint& start, const uint& end,
                                      const uint& chunk_size,
                                      const T& rates,
                                      pcg32& eng,
                                      std::exponential_distribution<double>& distr) {

    double r, key, X, w, t;

    // Create objects to store currently-selected position and key
    r = -1 * distr(eng);  // ~ log(U(0,1))
    key = r / rates[start]; // log(key)
    double largest_key = key;  // largest key (the one we're going to keep)
    uint largest_pos = start;  // position where largest key was found

    uint c = start;
    uint n_bp = rates.size(c);
    while (n_bp < chunk_size) {
        r = -1 * distr(eng);  // ~ log(U(0,1))
        X = r / largest_key;  // largest_key is already logged
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
                                const T& rates, pcg32& eng,
                                std::exponential_distribution<double>& distr) {

    double r, key, X, w, t;


    // Create objects to store currently-selected position and key
    r = -1 * distr(eng);        // ~ log(U(0,1))
    key = r / rates[start];     // log(key)
    double largest_key = key;   // largest key (the one we're going to keep)
    uint largest_pos = start;   // position where largest key was found

    uint c = start;
    while (c < end) {
        r = -1 * distr(eng);    // ~ log(U(0,1))
        X = r / largest_key;    // log(key)
        uint i = c + 1;
        X -= rates[c];
        w = rates[i];
        while (X > w && i < end) {
            X -= rates[i];
            i++;
            w = rates[i];
        }
        if (X > w) break;
        if (X <= 0) continue;

        largest_pos = i;

        t = std::exp(w * largest_key);  // key is log(key)
        r = runif_ab(eng, t, 1.0);
        key = std::log(r) / w;          // log(key)
        largest_key = key;

        c = i;

    }

    return largest_pos;
}



/*
 Simple class to wrap around (1) a vector of rates and (2) a vector of indices
 for a subsample of rates.
 Ultimately, this class allows you to iterate through indices for the indices vector
 using simply a bracket operator.
 For instance, if the rates vector is length 1e6 and we first subsample 1000 indices,
 we can do weighted sampling for just rates the 1000 indices refer to.
 This class is used to retrieve the rates based on indices from 0 to 999 instead of
 having to mess around with indices from 0 to 999999.
 */
template <typename T>
class RateGetter {
public:

    const T& rates;
    std::vector<uint> inds;

    RateGetter(const T& r, const uint& chunk_size)
        : rates(r), inds(chunk_size), distr(1.0) {};

    inline double operator[](const uint& idx) const {
        return rates[inds[idx]];
    }
    inline uint size() const noexcept {
        return inds.size();
    }

    inline void reset(pcg32& eng) {
        vitter_d<std::vector<uint>>(inds, rates.size(), eng);
    }

    inline double exp(pcg32& eng) {
        return distr(eng);
    }

private:
    std::exponential_distribution<double> distr;
};


template <typename T>
inline uint weighted_reservoir2_(RateGetter<T>& rates,
                                 pcg32& eng) {


    uint chunk = rates.size();

    rates.reset(eng);

    double r, key, X, w, t;


    // Create objects to store currently-selected position and key
    r = -1 * rates.exp(eng);    // ~ log(U(0,1))
    key = r / rates[0];         // log(key)
    double largest_key = key;   // largest key (the one we're going to keep)
    uint largest_pos = 0;       // position where largest key was found


    uint c = 0;
    while (c < (chunk-1)) {
        r = -1 * rates.exp(eng);    // ~ log(U(0,1))
        X = r / largest_key;        // log(key)
        uint i = c + 1;
        X -= rates[c];
        w = rates[i];
        while (X > w && i < (chunk-1)) {
            X -= rates[i];
            i++;
            w = rates[i];
        }
        if (X > w) break;
        if (X <= 0) continue;

        largest_pos = i;

        t = std::exp(w * largest_key);  // key is log(key)
        r = runif_ab(eng, t, 1.0);
        key = std::log(r) / w;          // log(key)
        largest_key = key;

        c = i;
    }

    return rates.inds[largest_pos];
}

#endif
