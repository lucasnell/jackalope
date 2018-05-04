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


//' C++ equivalent of R's \code{choose} function.
//'
//'
//' @param n Unsigned integer value. Make sure this isn't negative!
//' @param k Unsigned integer value. Make sure this isn't negative!
//'
//' @return Binomial coefficient (also integer).
//'
//' @noRd
//'
inline uint cpp_choose(const uint& n, uint k) {
    if (k > n) return 0;
    if (k * 2 > n) k = n - k;
    if (k == 0) return 1;

    uint result = n;
    for (uint i = 2; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}
// Same as above, but for doubles (even as doubles, they should be input as whole numbers)
inline double cpp_choose(const double& n, double k) {
    if (k > n) return 0;
    if (k * 2 > n) k = n - k;
    if (k == 0) return 1;

    double result = n;
    for (uint i = 2; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
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
