#ifndef __GEMINO_ALIAS_H
#define __GEMINO_ALIAS_H

/*
 ********************************************************

 Alias sampling

 ********************************************************
 */

#include <RcppArmadillo.h>
#include <sitmo.h>    // sitmo prng
#include <vector>  // vector class
#include <string>  // string class

#include "gemino_types.h" // integer types
#include "sequence_classes.h" // RefSequence class
#include "util.h"

using namespace Rcpp;


#define SMALL_TOLERANCE 0.00000001490116119385


namespace alias {
    const std::string nt_bases = "ACGT";
    const double sitmo_max = static_cast<double>(sitmo::prng_engine::max()) + 1;
}


// alias sampling of unsigned integers
class AliasUInts {
public:
    AliasUInts() : F(), L(), n(0) {};
    AliasUInts(const std::vector<double>& p, const double& tol = SMALL_TOLERANCE);
    // To get the length of F (and L bc they should always be the same)
    uint size() const noexcept {
        return n;
    }
    // Actual alias sampling
    inline uint sample(sitmo::prng_engine& eng) const {
        // uniform in range [0,1)
        double u = static_cast<double>(eng()) / alias::sitmo_max;
        // Not doing +1 [as is done in Yang (2006)] to keep it in 0-based indexing
        uint k = n * u;
        double r = n * u - k;
        if (r >= F[k]) k = L[k];
        return k;
    };
private:
    std::vector<double> F;
    std::vector<uint> L;
    uint n;
};





/*
 Alias sampling for indices in a `std::vector<uint>` or `arma::uvec`.
 It's assumed that the vector has been set to the desired size, and indexing
 is used to add values.
 */
template <typename T>
void alias_sample_uint(T& uints, const AliasUInts& sampler, sitmo::prng_engine& eng) {

    uint len = uints_get_size(uints);

    for (uint i = 0; i < len; i++) {
        uint k = sampler.sample(eng);
        uints[i] = k;
    }

    return;
}


/*
 Alias sampling for nucleotides in a string
 `nucleos` should already be set to the desired size using `nucleos = T(len, 'x')`
 or `nucleos.resize(len, 'x')`.
 `T` can be `std::string` or `RefSequence`.
 */
template <typename T>
void alias_sample_str(T& nucleos, const AliasUInts& sampler, sitmo::prng_engine& eng) {

    for (uint i = 0; i < nucleos.size(); i++) {
        uint k = sampler.sample(eng);
        nucleos[i] = alias::nt_bases[k];
    }

    return;
}


/*
 Class template for alias sampling a string, using an underlying AliasUInts object.
 `chars_in` should be the characters to sample from, `p` the probabilities of
 sampling those characters, and `tol` is the tolerance from the `AliasUInts` class
 construction.
 `T` can be `std::string` or `RefSequence`.
 */
template <typename T>
class AliasString {
public:

    T characters;

    AliasString(const T& chars_in, const std::vector<double>& p,
                const double& tol = SMALL_TOLERANCE)
        : characters(chars_in), uint_sampler(p, tol), n(p.size()) {
        if (p.size() != chars_in.size()) {
            stop("For an AliasString construction, arguments p and chars_in ",
                 "must be same length.");
        }
    }
    AliasString(const T& chars_in,
                const double& tol = SMALL_TOLERANCE)
        : characters(chars_in), uint_sampler(), n(chars_in.size()) {
        std::vector<double> p(n, 1 / static_cast<double>(n));
        uint_sampler(p, tol);
    }
    AliasString() {}

    std::string sample(const uint& N, sitmo::prng_engine& eng) const {
        std::string out(N, 'x');
        for (uint i = 0; i < N; i++) {
            uint k = uint_sampler.sample(eng);
            out[i] = characters[k];
        }
        return out;
    }

private:
    AliasUInts uint_sampler;
    uint n;
};


#endif
