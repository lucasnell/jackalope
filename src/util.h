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
#include <pcg/pcg_random.hpp> // pcg prng

#include "gemino_types.h"  // integer types
#include "pcg.h"  // runif_* methods


using namespace Rcpp;


/*
 Calling `base::options("width")$width`
 */
inline int get_width() {
    // Obtain environment containing function
    Environment base("package:base");
    // Make function callable from C++
    Function opt_r = base["options"];
    // Call the function and receive its list output
    List width_list = opt_r("width");
    int console_width = width_list["width"];
    return console_width;
}


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



// To return indices of sorted vector `values`.
template <typename T>
std::vector<uint32> increasing_indices(const std::vector<T>& values) {

    std::vector<uint32> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<uint32>(0));

    std::sort(
        begin(indices), end(indices),
        [&](size_t a, size_t b) { return values[a] < values[b]; }
    );
    return indices;
}
template <typename T>
std::vector<uint32> decreasing_indices(const std::vector<T>& values) {

    std::vector<uint32> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<uint32>(0));

    std::sort(
        begin(indices), end(indices),
        [&](size_t a, size_t b) { return values[a] > values[b]; }
    );
    return indices;
}




// Truncated normal when limit is not very far from the mean (< 5 SD away):
inline void trunc_rnorm_near(double& out,
                             const double& a_bar,
                             const double& mu,
                             const double& sigma,
                             pcg64& eng) {

    double p = R::pnorm5(a_bar, 0, 1, 1, 0);
    double u = runif_ab(eng, p, 1);

    double x = R::qnorm5(u, 0, 1, 1, 0);
    out = x * sigma + mu;

    return;
}
// And for when it IS very far from the mean:
inline void trunc_rnorm_far(double& out,
                            const double& a_bar,
                            const double& mu,
                            const double& sigma,
                            pcg64& eng) {
    double u, x_bar, v;
    u = runif_01(eng);
    x_bar = std::sqrt(a_bar * a_bar  - 2 * std::log(1 - u));
    v = runif_01(eng);
    while (v > (x_bar / a_bar)) {
        u = runif_01(eng);
        x_bar = std::sqrt(a_bar * a_bar  - 2 * std::log(1 - u));
        v = runif_01(eng);
    }
    out = sigma * x_bar + mu;

    return;
}




# endif
