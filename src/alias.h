#ifndef __GEMINO_ALIAS_H
#define __GEMINO_ALIAS_H

#include <RcppArmadillo.h>
#include <sitmo.h>    // sitmo prng
#include <vector>  // vector class
#include <string>  // string class
#include <unordered_map>  // unordered_map

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


// Alias sampling for indices
std::vector<uint> alias_sample_uint(const uint& N, const alias_FL& FL,
                                    sitmo::prng_engine& eng);

// Alias sampling for nucleotides in a string
void alias_sample_str(std::string& nucleos, const alias_FL& FL,
                      sitmo::prng_engine& eng);


#endif
