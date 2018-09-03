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
inline uint32 cpp_choose(const uint32& n, uint32 k) {
    if (k > n) return 0;
    if (k * 2 > n) k = n - k;
    if (k == 0) return 1;

    uint32 result = n;
    for (uint32 i = 2; i <= k; ++i) {
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
    for (uint32 i = 2; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

/*
 Get a size from either an arma::uvec or std::vector<uint32>.
 This is used in template functions that work for either class.
 */
inline uint32 uints_get_size(std::vector<uint32>& uints) {
    return uints.size();
}
inline uint32 uints_get_size(arma::uvec& uints) {
    return uints.n_elem;
}


//' GC proportion of a single string.
//'
//'
//' @param sequence String for a single sequence.
//'
//' @return Proportion of sequence that's a `'G'` or `'C'`.
//'
//' @noRd
//'
inline double gc_prop(const std::string& sequence) {
    double total_seq = sequence.size();
    double total_gc = 0;
    for (uint32 i = 0; i < total_seq; i++) {
        if (sequence[i] == 'G' || sequence[i] == 'C') {
            total_gc += 1;
        }
    }
    double gc_prop = total_gc / total_seq;
    return gc_prop;
}
// ... overloaded for portion of a string
inline double gc_prop(const std::string& sequence,
                      const uint32& start,
                      const uint32& stop) {
    double total_seq = stop - start + 1;
    double total_gc = 0;
    for (uint32 i = start; i <= stop; i++) {
        if (sequence[i] == 'G' || sequence[i] == 'C') {
            total_gc += 1;
        }
    }
    double gc_prop = total_gc / total_seq;
    return gc_prop;
}



# endif
