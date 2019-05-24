# ifndef __JACKAL_UTIL_H
# define __JACKAL_UTIL_H


/*
 ********************************************************

 Miscellaneous helper functions.

 ********************************************************
 */



#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <pcg/pcg_random.hpp> // pcg prng

#ifdef _OPENMP
#include <omp.h>  // omp
#endif


#include "jackalope_types.h"  // integer types
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



// Makes integers print with commas
// Do NOT use this with doubles
template <typename T>
inline std::string big_int_format(const T& x) {
    std::string y = std::to_string(x);
    if (y.size() > 3) {
        uint64 i = 3;
        while (i < y.size()) {
            y.insert(y.end()-i, ',');
            i += 3 + 1;
        }
    }
    return y;
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



//' GC proportion for a single string.
//'
//'
//' @param sequence String for a single sequence.
//' @param start Position in sequence at which to begin the calculation.
//' @param stop Position in sequence at which to end the calculation.
//'
//' @return Proportion of sequence that's a `'G'` or `'C'`.
//'
//' @noRd
//'
inline double gc_prop(const std::string& sequence) {
    double total_seq = sequence.size();
    double total_gc = 0;
    for (uint64 i = 0; i < total_seq; i++) {
        if (sequence[i] == 'G' || sequence[i] == 'C') {
            total_gc += 1;
        }
    }
    double gc_prop = total_gc / total_seq;
    return gc_prop;
}
inline double gc_prop(const std::string& sequence,
                      const uint64& start,
                      const uint64& stop) {
    double total_seq = stop - start + 1;
    double total_gc = 0;
    for (uint64 i = start; i <= stop; i++) {
        if (sequence[i] == 'G' || sequence[i] == 'C') {
            total_gc += 1;
        }
    }
    double gc_prop = total_gc / total_seq;
    return gc_prop;
}


//' Same as above, but for any one nucleotide.
//'
//'
//' @inheritParams gc_prop
//' @param nt Character to count for the proportion.
//'
//' @return Proportion of sequence that's a `nt`.
//'
//' @noRd
//'
inline double nt_prop(const std::string& sequence,
                      const char& nt) {
    double total_seq = sequence.size();
    double total_nt = 0;
    for (uint64 i = 0; i < total_seq; i++) {
        if (sequence[i] == nt) total_nt += 1;
    }
    double nt_prop = total_nt / total_seq;
    return nt_prop;
}
// Overloaded for part of a sequence
inline double nt_prop(const std::string& sequence,
                      const char& nt,
                      const uint64& start,
                      const uint64& stop) {
    double total_seq = stop - start + 1;
    double total_nt = 0;
    for (uint64 i = start; i <= stop; i++) {
        if (sequence[i] == nt) total_nt += 1;
    }
    double nt_prop = total_nt / total_seq;
    return nt_prop;
}





//' For prettier long error messages.
//'
//'
inline void str_stop(const std::vector<std::string>& err_msg_vec) {
    std::string err_msg = "";
    for (const std::string& err : err_msg_vec) err_msg += err;
    throw(Rcpp::exception(err_msg.c_str(), false));
}
//' For prettier long warning messages.
//'
//'
inline void str_warn(const std::vector<std::string>& warn_msg_vec) {
    std::string warn_msg = "";
    for (const std::string& warn : warn_msg_vec) warn_msg += warn;
    Rcpp::warning(warn_msg.c_str());
    return;
}






//' Check that the number of threads doesn't exceed the number available, and change
//' to 1 if OpenMP isn't enabled.
//'
//' @noRd
//'
inline void thread_check(uint64& n_threads) {

#ifdef _OPENMP
    if (n_threads == 0) n_threads = 1;
    if (n_threads > omp_get_max_threads()) {
        std::string max_threads = std::to_string(omp_get_max_threads());
        str_stop({"\nThe number of requested threads (", std::to_string(n_threads),
                 ") exceeds the max available on the system (", max_threads, ")."});
    }
#else
    n_threads = 1;
#endif

    return;
}





//' Split integer `x` into `n` chunks that are as even as possible.
//'
//' @noRd
//'
inline std::vector<uint64> split_int(const uint64& x,
                                     const uint64& n) {

    std::vector<uint64> out(n, x / n);
    uint64 sum_reads = n * static_cast<uint64>(x / n);
    uint64 i = 0;
    while (sum_reads < x) {
        out[i]++;
        i++;
        sum_reads++;
    }
    return out;

}





# endif
