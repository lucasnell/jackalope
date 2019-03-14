#include <RcppArmadillo.h>
#include <vector>
#include <string>

#include "gemino_types.hpp"  // integer types


using namespace Rcpp;

//[[Rcpp::export]]
bool using_openmp() {
    bool out = false;
#ifdef _OPENMP
    out = true;
#endif
    return out;
}


