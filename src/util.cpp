


#include <RcppArmadillo.h>
#include <vector>
#include <string>

#include "gemino_types.h"  // integer types
#include "util.h"

using namespace Rcpp;








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
    int total_seq = sequence.size();
    double total_gc = 0;
    for (int i = 0; i < total_seq; i++) {
        if (sequence[i] == 'G' || sequence[i] == 'C') {
            total_gc += 1;
        }
    }
    double gc_prop = total_gc / total_seq;
    return gc_prop;
}
// ... overloaded for portion of a string
double gc_prop(const std::string& sequence, const uint& start, const uint& stop) {
    int total_seq = stop - start + 1;
    double total_gc = 0;
    for (int i = start; i <= stop; i++) {
        if (sequence[i] == 'G' || sequence[i] == 'C') {
            total_gc += 1;
        }
    }
    double gc_prop = total_gc / total_seq;
    return gc_prop;
}
