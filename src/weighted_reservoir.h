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
 Templates that do weighted reservoir sampling for a vector, returning
 one unsigned integer index to the sampled location in the `rates` vector.

 The first version simply samples along the entire vector.

 The second, "chunk" version first sub-samples (*not* weighted) a set number of
 indices, then only does weighted sampling on rates with those indices.
 This can improve performance pretty significantly, depending on the number of
 sub-samples compared to the whole vector's size.

 Both versions require that class `T` has...
 (1) a bracket-operator method (`operator[]`) that returns the rate (type `double`)
     for a given index
 (2) a `size()` method that returns the max index to sample + 1.

 Both versions also require that you pre-construct an object to retrieve the rates
 from. The classes you need to construct are different for the first and second
 sampling functions below, and they are templates with the same requirements for
 class `T` as listed above.
 See the classes `ReservoirRates` and `ChunkReservoirRates` below for more info.

 */




/*
 Regular weighted reservoir sampling.
 This template does most of the work for the two `sample` methods for the two classes
 in this file.
 */

template <typename R>
inline uint weighted_reservoir_(R& obj, pcg32& eng) {

    double r, key, X, w, t;

    uint end = obj.rates.size() - 1;

    // Create objects to store currently-selected position and key
    r = -1 * obj.rexp_(eng);        // ~ log(U(0,1))
    key = r / obj.rates[0];     // log(key)
    double largest_key = key;   // largest key (the one we're going to keep)
    uint largest_pos = 0;   // position where largest key was found

    uint c = 0;
    while (c < end) {
        r = -1 * obj.rexp_(eng);    // ~ log(U(0,1))
        X = r / largest_key;        // log(key)
        uint i = c + 1;
        X -= obj.rates[c];
        w = obj.rates[i];
        while (X > w && i < end) {
            X -= obj.rates[i];
            i++;
            w = obj.rates[i];
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
 A class for weighted reservoir sampling.
 This is just a wrapper for a reference to a vector of rates plus an exponential
 distribution object.
 I made this class to make sure the exponential distribution isn't messed with, plus
 to be consistent with the "chunked" version below.
 Class `T` needs to return a double with square brackets and have a `size()` method.
 */
template <typename T>
class ReservoirRates {

public:

    T rates;

    ReservoirRates() : rates(), distr(1.0) {};
    ReservoirRates(const T& r) : rates(r), distr(1.0) {};
    /*
     `cs` is ignored here, and this constructor is only to allow for template use along
     with chunked version:
     */
    ReservoirRates(const T& r, const uint& cs) : rates(r), distr(1.0) {};
    // Assignment operator
    ReservoirRates& operator=(const ReservoirRates& other) {
        rates = other.rates;
    }

    inline double rexp_(pcg32& eng) {
        return distr(eng);
    }

    // Sample for one location
    inline uint sample(pcg32& eng) {
        return weighted_reservoir_<ReservoirRates, T>(*this, eng);
    }

private:
    std::exponential_distribution<double> distr;

};







/*
 A class for "chunked" weighted reservoir sampling.
 Simple class to wrap around (1) a vector of rates and (2) a vector of indices
 for a subsample of rates. (The former is a reference to prevent unnecessary copying.)
 Ultimately, this class allows you to iterate through indices for the indices vector
 using simply a bracket operator.
 For instance, if the rates vector is length 1e6 and we first subsample 1000 indices,
 we can do weighted sampling for just rates the 1000 indices refer to.
 This class is used to retrieve the rates based on indices from 0 to 999 instead of
 having to mess around with indices from 0 to 999999.

 Class `T` needs to return a double with square brackets and have a `size()` method.
 */

template <typename T>
struct ChunkRateGetter {

    T all_rates;
    std::vector<uint> inds;

    ChunkRateGetter() : all_rates(), inds() {};
    ChunkRateGetter(const T& r, const uint& chunk)
        : all_rates(r), inds(chunk) {};
    // Assignment operator
    ChunkRateGetter& operator=(const ChunkRateGetter& other) {
        all_rates = other.all_rates;
        inds = other.inds;
    }

    inline double operator[](const uint& idx) const {
        return all_rates[inds[idx]];
    }
    inline uint size() const noexcept {
        return inds.size();
    }
    inline void reset(pcg32& eng) {
        vitter_d<std::vector<uint>>(inds, all_rates.size(), eng);
    }

    /*
     Return index referring to a position in `all_rates` using an index for `inds`.
     For example, if `all_rates` is length 200 and `inds` is length 100, the input
     index here should never be > 99, but the output index could be up to 199.
     */
    inline uint operator()(const uint& idx) const {
        return inds[idx];
    }
};

template <typename T>
class ChunkReservoirRates {

public:

    ChunkRateGetter<T> rates;

    ChunkReservoirRates() : rates(), distr(1.0) {};
    ChunkReservoirRates(const T& r, const uint& chunk)
        : rates(r, chunk), distr(1.0) {};
    // Assignment operator
    ChunkReservoirRates<T>& operator=(const ChunkReservoirRates<T>& other) {
        rates = other.rates;
        distr = std::exponential_distribution<double>(1.0);
    }

    inline double rexp_(pcg32& eng) {
        return distr(eng);
    }

    // Sample for one location
    inline uint sample(pcg32& eng) {
        rates.reset(eng);
        uint i = weighted_reservoir_<ChunkReservoirRates, T>(*this, eng);
        return rates(i);
    }


private:

    std::exponential_distribution<double> distr;

};




#endif
