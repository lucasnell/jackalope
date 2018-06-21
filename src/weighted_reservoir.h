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
 (2) a `size()` method that returns the max index to sample + 1

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

 If `end` is set to zero and `start >= end`, then `end` will be reset to the maximum
 possible value.
 This happens by default.
 */

template <typename R>
inline uint32 weighted_reservoir_(R& obj, pcg32& eng,
                                const uint32& start = 0,
                                uint32 end = 0) {

    double r, key, X, w, t;

    if (start >= end && end == 0) end = obj.res_rates.size() - 1;

    // Create objects to store currently-selected position and key
    r = -1 * obj.rexp_(eng);        // ~ log(U(0,1))
    key = r / obj.res_rates[start];     // log(key)
    double largest_key = key;   // largest key (the one we're going to keep)
    uint32 largest_pos = start;   // position where largest key was found

    uint32 c = start;
    while (c < end) {
        r = -1 * obj.rexp_(eng);    // ~ log(U(0,1))
        X = r / largest_key;        // log(key)
        uint32 i = c + 1;
        X -= obj.res_rates[c];
        w = obj.res_rates[i];
        while (X > w && i < end) {
            X -= obj.res_rates[i];
            i++;
            w = obj.res_rates[i];
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

    T res_rates;

    ReservoirRates() : res_rates(), distr(1.0) {};
    ReservoirRates(const T& r) : res_rates(r), distr(1.0) {};
    ReservoirRates(const ReservoirRates<T>& other)
        : res_rates(other.res_rates), distr(1.0) {}
    /*
     `cs` is ignored here, and this constructor is only to allow for template use along
     with chunked version:
     */
    ReservoirRates(const T& r, const uint32& cs) : res_rates(r), distr(1.0) {};
    // Assignment operator
    ReservoirRates<T>& operator=(const ReservoirRates<T>& other) {
        res_rates = other.res_rates;
        return *this;
    }

    inline double rexp_(pcg32& eng) {
        return distr(eng);
    }

    // Sample for one location
    inline uint32 sample(pcg32& eng) {
        return weighted_reservoir_<ReservoirRates<T>>(*this, eng);
    }
    // Sample for one location inside a range
    inline uint32 sample(pcg32& eng, const uint32& start, const uint32& end) {
        return weighted_reservoir_<ReservoirRates<T>>(*this, eng, start, end);
    }

protected:
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
class ChunkRateGetter {

public:

    // The vector of rates
    T all_rates;
    // The vector of sub-sampled indices (possible values from 0 to all_rates.size() - 1)
    std::vector<uint32> inds;

    // To store the initial chunk size argument to the constructor.
    uint32 chunk_size;


    ChunkRateGetter() : all_rates(), inds(), chunk_size(0), use_vitter(true) {};
    ChunkRateGetter(const T& r, const uint32& chunk_size_)
        : all_rates(r), inds(chunk_size_), chunk_size(chunk_size_), use_vitter(true) {
        /*
         If chunk >= all_rates size, then turn `inds` into a vector from 0 to
         `all_rates.size() - 1`, and set use_vitter to false.
         This makes it easier to avoid the uniform sampling step when using
         this class.
         */
        if (chunk_size >= all_rates.size()) {
            inds.resize(all_rates.size());
            for (uint32 i = 0; i < all_rates.size(); i++) inds[i] = i;
            use_vitter = false;
        }
    };
    ChunkRateGetter(const ChunkRateGetter<T>& other)
        : all_rates(other.all_rates), inds(other.inds), chunk_size(other.chunk_size),
          use_vitter(other.use_vitter) {}
    // Assignment operator
    ChunkRateGetter<T>& operator=(const ChunkRateGetter<T>& other) {
        all_rates = other.all_rates;
        inds = other.inds;
        chunk_size = other.chunk_size;
        use_vitter = other.use_vitter;
        return *this;
    }

    inline double operator[](const uint32& idx) const {
        return all_rates[inds[idx]];
    }
    inline uint32 size() const noexcept {
        return inds.size();
    }
    void reset(pcg32& eng) {
        recheck_size_();
        // Skip vitter_d if `inds` is a vector from 0 to `all_rates.size() - 1`:
        if (!use_vitter) return;
        // Otherwise, sample uniformly:
        vitter_d<std::vector<uint32>>(inds, all_rates.size(), eng);
        return;
    }
    void reset(pcg32& eng, const uint32& start, const uint32& end) {

        if (start >= all_rates.size() || start > end) {
            stop("start too high in ChunkRateGetter::reset");
        }
        if (end >= all_rates.size()) stop("end too high in ChunkRateGetter::reset");

        uint32 range_size = end - start + 1;

        if (range_size <= chunk_size) {
            adjust_range_small_(start, end, range_size);
        } else {
            adjust_range_big_(start, end, range_size);
        }

        // Skip vitter_d if `inds` is a vector from `start` to `end`:
        if (!use_vitter) return;
        // Otherwise, sample uniformly:
        vitter_d<std::vector<uint32>>(inds, range_size, eng, start);
        return;
    }



private:


    // Boolean for whether to use vitter_d first to do uniform sampling.
    bool use_vitter;


    /*
     Check to see if...
     1) something (e.g., a deletion) has changed `all_rates` to make it go from larger
     to smaller than `chunk_size`
     2) something (e.g., an insertion) has changed `all_rates` to make it go from
     smaller to larger than `inds`
     Use this for unranged sampling only!
     */
    void recheck_size_() {
        /*
         If `all_rates.size()` is now <= `chunk_size` but `use_vitter` is still true,
         then make `inds` a vector from 0 to `all_rates.size() - 1`
         and set `use_vitter` to false
         */
        if (all_rates.size() <= chunk_size) {
            if (use_vitter) {
                inds.resize(all_rates.size());
                for (uint32 i = 0; i < all_rates.size(); i++) inds[i] = i;
                use_vitter = false;
            }
            return;
        }
        /*
         If `all_rates.size()` is now > `chunk_size` but `use_vitter` is still false,
         then make `inds` a vector of length `chunk_size`
         and set `use_vitter` to true
         */
        if (!use_vitter) {
            inds.resize(chunk_size);
            use_vitter = true;
        }

        return;
    }

    /*
     These two functions are for if you want to adjust this object for ranged
     sampling when the range endpoints and `all_rates` change.
     It's important that you explicitly change the range, otherwise you may cause
     weird errors when doing a ranged reset.
     Use these for ranged sampling only!
     */

    // For when range_size <= chunk_size
    inline void adjust_range_small_(const uint32& start, const uint32& end,
                                    const uint32& range_size) {

        /*
         Note #1:
         Because `range_size < chunk_size`, the size of the output will be `range_size`.

         Note #2:
         Because we are using a range that should reference `all_rates`,
         `range_size` should never be > `all_rates.size()`, but it can be equal to it.
         There are already checks to ensure `range_size <= all_rates.size()`,
         so I won't put an extra one here.

         Note #3:
         Because `range_size <= chunk_size` (based on a check done before running this
         function) and because we only want to sample within the range,
         we don't want to sub-sample first.
         So if we run this function only when `range_size <= chunk_size` (as we should!),
         `use_vitter` should always be false and `inds` should go from `start` to `end`.
         The rest of this function is to try to efficiently manipulate `inds` properly.
         */

        uint32 inds_size = inds.size();

        /*
         If `use_vitter` was false, then `inds` should currently be a vector from `x` to
         `inds.size() + x - 1`, where `x` was a previous `start` value.
         I just need to adjust this so it now goes from `start` to `end`.
         I'm assuming that if `x == start`, then the whole vector `inds` is the same
         as the range `start` to `end` if `inds.size() == range_size`.
        */
        if (!use_vitter) {
            if (inds_size < range_size) {
                inds.reserve(range_size);
                if (inds_size == 0) {
                    for (uint32 i = start; i <= end; i++) inds.push_back(i);
                    return;
                }
                if (inds.front() != start) {
                    for (uint32 i = 0; i < inds_size; i++) inds[i] = start + i;
                }
                for (uint32 i = inds_size + start; i <= end; i++) inds.push_back(i);
            } else {
                if (inds_size > range_size) inds.resize(range_size);
                if (inds.front() != start) {
                    for (uint32 i = 0; i < range_size; i++) inds[i] = start + i;
                }
            }
            /*
             If `use_vitter` was true, then `inds` should NOT be a vector from
            `x` to `inds.size() + x - 1`, where `x` was a previous `start` value
            */
        } else {
            use_vitter = false;
            if (inds_size < range_size) {
                inds.reserve(range_size);
                // (The check for `inds_size == 0` is not needed here.)
                for (uint32 i = 0; i < inds_size; i++) inds[i] = start + i;
                for (uint32 i = inds_size + start; i <= end; i++) inds.push_back(i);
            } else {
                if (inds_size > range_size) inds.resize(range_size);
                for (uint32 i = 0; i < range_size; i++) inds[i] = start + i;
            }
        }

        return;
    }


    // For when range_size > chunk_size
    inline void adjust_range_big_(const uint32& start, const uint32& end,
                                  const uint32& range_size) {

        /*
         Note #1:
         Because `range_size > chunk_size`, the size of the output will be `chunk_size`.

         Note #2:
         Because we've already made sure that `range_size <= all_rates.size()`
         and `range_size > chunk_size`,
         we also know that `chunk_size < all_rates.size()`.
         For this reason, `use_vitter` should always be set to true when this function
         is run.
         */

        if (inds.size() != chunk_size) inds.resize(chunk_size);

        use_vitter = true;

        return;
    }

};


template <typename T>
class ChunkReservoirRates {

public:

    ChunkRateGetter<T> res_rates;

    ChunkReservoirRates() : res_rates(), distr(1.0) {};
    ChunkReservoirRates(const T& r, const uint32& chunk)
        : res_rates(r, chunk), distr(1.0) {};
    ChunkReservoirRates(const ChunkReservoirRates<T>& other)
        : res_rates(other.res_rates), distr(1.0) {}
    // Assignment operator
    ChunkReservoirRates<T>& operator=(const ChunkReservoirRates<T>& other) {
        res_rates = other.res_rates;
        return *this;
    }

    inline double rexp_(pcg32& eng) {
        return distr(eng);
    }

    // Sample for one location
    inline uint32 sample(pcg32& eng) {
        res_rates.reset(eng);
        uint32 i = weighted_reservoir_<ChunkReservoirRates<T>>(*this, eng);
        return res_rates.inds[i];
    }
    // Sample for one location inside a range
    inline uint32 sample(pcg32& eng, const uint32& start, const uint32& end) {
        res_rates.reset(eng, start, end);
        uint32 i = weighted_reservoir_<ChunkReservoirRates<T>>(*this, eng);
        return res_rates.inds[i];
    }

protected:

    std::exponential_distribution<double> distr;

};




#endif
