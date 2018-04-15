


#include <RcppArmadillo.h>
#include <vector>
#include <string>

#include "gemino_types.h"  // integer types
#include "util.h"

using namespace Rcpp;



/*
 Clear memory from a std::vector, std::deque, or std::string.
 Simply erasing objects does not clear memory.
 I'm not sure this works for other classes. You should check before using this function.
 */
template <class U>
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
