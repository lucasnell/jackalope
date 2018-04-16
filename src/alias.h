#ifndef __GEMINO_ALIAS_H
#define __GEMINO_ALIAS_H

#include <RcppArmadillo.h>
#include <sitmo.h>    // sitmo prng
#include <vector>  // vector class
#include <string>  // string class

#include "gemino_types.h" // integer types
#include "sequence_classes.h" // RefSequence class
#include "util.h"

using namespace Rcpp;


namespace alias {
    const std::string nt_bases = "ACGT";
    const double sitmo_max = static_cast<double>(sitmo::prng_engine::max()) + 1;
}


// F and L vectors for alias sampling
class alias_FL {
public:
    std::vector<double> F;
    std::vector<uint> L;
    alias_FL() : F(), L() {};
    alias_FL(const std::vector<double>& p, const double& tol = 0.00000001490116119385);
    // To get the length of F (and L bc they should always be the same)
    uint size() const noexcept {
        return F.size();
    }
};



// Actual alias sampling
inline uint alias_sample(const uint& n, const alias_FL& FL,
                         sitmo::prng_engine& eng) {
    // uniform in range [0,1)
    double u = static_cast<double>(eng()) / alias::sitmo_max;
    // Not doing +1 [as is done in Yang (2006)] to keep it in 0-based indexing
    uint k = n * u;
    double r = n * u - k;
    if (r >= FL.F[k]) k = FL.L[k];
    return k;
}





/*
 Alias sampling for indices in a `std::vector<uint>` or `arma::uvec`.
 It's assumed that the vector has been set to the desired size, and indexing
 is used to add values.
 */
template <typename T>
void alias_sample_uint(T& uints, const alias_FL& FL, sitmo::prng_engine& eng) {

    uint n = FL.size();

    uint len = uints_get_size(uints);

    for (uint i = 0; i < len; i++) {
        uint k = alias_sample(n, FL, eng);
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
void alias_sample_str(T& nucleos, const alias_FL& FL, sitmo::prng_engine& eng) {

    uint n = FL.size();
    if (n != 4) stop("Function alias_sample_str only uses alias_FL objects of size 4.");

    for (uint i = 0; i < nucleos.size(); i++) {
        uint k = alias_sample(n, FL, eng);
        nucleos[i] = alias::nt_bases[k];
    }

    return;
}




#endif
