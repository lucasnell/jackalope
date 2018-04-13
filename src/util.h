# ifndef __GEMINO_UTIL_H
# define __GEMINO_UTIL_H

#include <RcppArmadillo.h>
#include <vector>
#include <string>

#include "gemino_types.h"  // integer types


using namespace Rcpp;



std::string cpp_rando_seq(uint len);


double gc_prop(const std::string& sequence);
double gc_prop(const std::string& sequence, const uint& start, const uint& stop);

# endif
