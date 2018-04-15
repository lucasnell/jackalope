#ifndef __GEMINO_ALIAS_H
#define __GEMINO_ALIAS_H

#include <RcppArmadillo.h>
#include <sitmo.h>    // sitmo prng
#include <vector>  // vector class
#include <string>  // string class

#include "gemino_types.h" // integer types


using namespace Rcpp;


namespace alias {
    const std::string nt_bases = "ACGT";
    const double alias_sitmo_max = static_cast<double>(sitmo::prng_engine::max()) + 1;
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


/*
 Alias sampling for indices.
 `U` can be `std::vector<uint>` or `arma::uvec`
 */
template <class U>
void alias_sample_uint(U& uints, const alias_FL& FL, sitmo::prng_engine& eng);

/*
 Alias sampling for nucleotides in a string
 `nucleos` should already be set to the desired size using `nucleos = U(len, 'x')`
 or `nucleos.resize(len, 'x')`.
 `U` can be `std::string` or `RefSequence`.
 */
template <class U>
void alias_sample_str(U& nucleos, const alias_FL& FL, sitmo::prng_engine& eng);


#endif
