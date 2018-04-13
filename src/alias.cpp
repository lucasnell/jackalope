
/*
 ========================================================================================
 ========================================================================================

 Alias sampling: much more efficient way to do weighted sampling.

 For more info, see p 299 in...

 Yang, Z. 2006. Computational molecular evolution. (P. H. Harvey and R. M. May, Eds.).
 Oxford University Press, New York, NY, USA.

 ========================================================================================
 ========================================================================================
 */



#include <RcppArmadillo.h>
#include <sitmo.h>    // sitmo prng
#include <vector>  // vector class
#include <string>  // string class
#include <unordered_map>  // unordered_map

#include "gemino_types.h" // integer types
#include "alias.h" // alias sampling
#include "util.h"

using namespace Rcpp;



// Construct the L and F vectors in a alias_FL object

alias_FL::alias_FL(const std::vector<double>& probs, const double& tol) {

    uint n = probs.size();

    // F_ and L_ are temporary vectors in arma formats to make the math easier
    arma::vec p(probs);
    double sum_p = arma::accu(p);
    p /= sum_p;
    arma::vec F_ = n * p;
    arma::uvec L_ = arma::regspace<arma::uvec>(0, n - 1);
    arma::ivec I(n);
    for (uint i = 0; i < n; i++) {
        if (F_(i) == 1) {
            L_(i) = i;
            I(i) = 0;
        } else if (F_(i) < 1) {
            I(i) = -1;
        } else {
            I(i) = 1;
        }
    }

    while (arma::any(I != 0)) {

        arma::uvec jv = arma::find(I == -1);  // underfull (i.e., F_ < 1)
        arma::uvec kv = arma::find(I == 1);  // overfull (i.e. F_ > 1)
        uint j = jv(0);
        if (kv.n_elem == 0) {
            stop("Numerical issue. Difference between one of the entries ",
                 "and 1 is " + std::to_string(F_(j) - 1));
        }
        uint k = kv(0);
        L_(j) = k;
        F_(k) = F_(k) - (1 - F_(j));
        I(j) = 0;
        if (std::abs(1 - F_(k)) < tol) {
            F_(k) = 1;
            I(k) = 0;
        } else if (F_(k) < 1) {
            I(k) = -1;
        }
    }

    F = arma::conv_to<std::vector<double>>::from(F_);
    L = arma::conv_to<std::vector<uint>>::from(L_);

    return;

}



// Alias sampling for indices
std::vector<uint> alias_sample_uint(const uint& N,
                                    const alias_FL& FL,
                                    sitmo::prng_engine& eng) {

    std::vector<uint> out(N);

    uint n = FL.size();

    for (uint i = 0; i < N; i++) {
        // uniform in range [0,1)
        double u = static_cast<double>(eng()) / alias::alias_sitmo_max;
        // Not doing +1 [as is done in Yang (2006)] to keep it in 0-based indexing
        uint k = n * u;
        double r = n * u - k;
        if (r >= FL.F[k]) k = FL.L[k];
        out[i] = k;
    }

    return out;
}

// Alias sampling for nucleotides in a string
// `nucleos` should already be allocated for the desired size
void alias_sample_str(std::string& nucleos,
                      const alias_FL& FL,
                      sitmo::prng_engine& eng) {

    uint N = nucleos.size();

    uint n = FL.size();
    if (n != 4) stop("Function alias_sample_str only uses alias_FL objects of size 4.");

    for (uint i = 0; i < N; i++) {
        // uniform in range [0,1)
        double u = static_cast<double>(eng()) / alias::alias_sitmo_max;
        // Not doing +1 [as is done in Yang (2006)] to keep it in 0-based indexing
        uint k = n * u;
        double r = n * u - k;
        if (r >= FL.F[k]) k = FL.L[k];
        nucleos[i] = alias::nt_bases[k];
    }

    return;
}

