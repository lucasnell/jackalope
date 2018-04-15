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


/*
 Clear memory from a std::vector, std::deque, or std::string.
 Simply erasing objects does not clear memory.
 I'm not sure this works for other classes. You should check before using this function.
 */
template <typename U>
void clear_memory(U& x) {
    U(x.begin(), x.end()).swap(x);
}


/*
 Get a size from either an arma::uvec or std::vector<uint>.
 This is used in template functions that work for either class.
 */
inline uint uints_get_size(std::vector<uint>& uints) {
    return uints.size();
}
inline uint uints_get_size(arma::uvec& uints) {
    return uints.n_elem;
}


double gc_prop(const std::string& sequence);
double gc_prop(const std::string& sequence, const uint& start, const uint& stop);




# endif
