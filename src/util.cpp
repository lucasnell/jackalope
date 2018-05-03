


#include <RcppArmadillo.h>
#include <vector>
#include <string>

#include "gemino_types.h"  // integer types
#include "util.h"

using namespace Rcpp;


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
    if (k > n) {
        return 0;
    }
    if (k * 2 > n) {
        k = n - k;
    }
    if (k == 0) {
        return 1;
    }

    uint result = n;
    for (uint i = 2; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
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
double gc_prop(const std::string& sequence) {
    double total_seq = sequence.size();
    double total_gc = 0;
    for (uint i = 0; i < total_seq; i++) {
        if (sequence[i] == 'G' || sequence[i] == 'C') {
            total_gc += 1;
        }
    }
    double gc_prop = total_gc / total_seq;
    return gc_prop;
}
// ... overloaded for portion of a string
double gc_prop(const std::string& sequence, const uint& start, const uint& stop) {
    double total_seq = stop - start + 1;
    double total_gc = 0;
    for (uint i = start; i <= stop; i++) {
        if (sequence[i] == 'G' || sequence[i] == 'C') {
            total_gc += 1;
        }
    }
    double gc_prop = total_gc / total_seq;
    return gc_prop;
}
