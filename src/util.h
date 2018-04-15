# ifndef __GEMINO_UTIL_H
# define __GEMINO_UTIL_H


/*
 ********************************************************

 Miscellaneous helper functions.

 ********************************************************
 */



#include <RcppArmadillo.h>
#include <vector>
#include <string>

#include "gemino_types.h"  // integer types


using namespace Rcpp;


template <class U>
void clear_memory(U& x);

inline uint uints_get_size(std::vector<uint>& uints);
inline uint uints_get_size(arma::uvec& uints);

double gc_prop(const std::string& sequence);
double gc_prop(const std::string& sequence, const uint& start, const uint& stop);




# endif
