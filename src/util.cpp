#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <vector>
#include <string>

#include "jackalope_types.h"  // integer types


using namespace Rcpp;

//[[Rcpp::export]]
bool using_openmp() {
    bool out = false;
#ifdef _OPENMP
    out = true;
#endif
    return out;
}


