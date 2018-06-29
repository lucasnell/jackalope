
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

#include "gemino_types.h" // for integer types
#include "str_manip.h" // for lookup tables

using namespace Rcpp;


/*
 Filter for only T, C, A, G, N, t, c, a, g, or n. Others characters are ignored.
 If upper=true, it converts lowercase to uppercase.
 */
void filter_nucleos(std::string& nucleos, const bool& upper) {
    if (upper) {
        for (char& c : nucleos) c = str_manip::upper_filter_table[c];
    } else {
        for (char& c : nucleos) c = str_manip::filter_table[c];
    }
    return;
}



/*
 Split a string based on a single-character delimiter
 */

std::vector<std::string> cpp_str_split_delim(const std::string& in_string,
                                             const char& split) {


    std::vector<std::string> out(1, "");
    std::string::size_type n = 1;

    std::string::size_type i = in_string.find(split);
    if (i != std::string::npos) {
        // Index for the output vector
        uint32 j = 0;
        // Index for the previous i:
        uint32 i0 = 0;
        while (i != std::string::npos) {
            for (std::string::size_type k = i0; k < i; k++) {
                out[j] += in_string[k];
            }
            i0 = i + n;
            i = in_string.find(split, i0);
            j++;
            out.push_back("");
        }
        for (std::string::size_type k = i0; k < in_string.size(); k++) {
            out[j] += in_string[k];
        }
    } else {
        out[0] = in_string;
    }

    return out;
}
