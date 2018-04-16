#include <RcppArmadillo.h>
#include <sitmo.h>    // sitmo prng
#include <vector>  // vector class
#include <string>  // string class

#include "gemino_types.h" // integer types
#include "sequence_classes.h" // RefSequence class
#include "util.h"
#include "alias.h"

using namespace Rcpp;



// Construct the L and F vectors in a AliasUInts object

AliasUInts::AliasUInts(const std::vector<double>& probs, const double& tol) {

    n = probs.size();

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



