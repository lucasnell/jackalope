# ifndef __JACKALOPE_STR_MANIP_H
# define __JACKALOPE_STR_MANIP_H


#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

#include "jackalope_types.h" // integer types

using namespace Rcpp;

namespace str_manip {

/*
 Lookup table to turn lowercase to uppercase, but only for T, C, A, G, and N (plus
 lowercase versions). Anything else is a zero, which will help make sure
 nothing weird gets read.
 */
const std::vector<uint64> upper_filter_table = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 65, 0, 67, 0, 0, 0, 71, 0, 0, 0, 0, 0, 0, 78, 0,
    0, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 65, 0, 67,
    0, 0, 0, 71, 0, 0, 0, 0, 0, 0, 78, 0, 0, 0, 0, 0, 84, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

// Same thing, but just filtering, no converting to uppercase
const std::vector<uint64> filter_table = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 65, 0, 67, 0, 0, 0, 71, 0, 0, 0, 0, 0, 0, 78, 0,
    0, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 97, 0, 99,
    0, 0, 0, 103, 0, 0, 0, 0, 0, 0, 110, 0, 0, 0, 0, 0, 116, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

// For complements (Ns stay as Ns)
const std::vector<uint64> cmp_map = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 84, 0, 71, 0, 0, 0, 67, 0, 0, 0, 0, 0, 0, 78, 0,
    0, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

}


/*
 Filter for only T, C, A, G, N, t, c, a, g, or n. Others characters are ignored.
 If upper=true, it converts lowercase to uppercase.
 */
inline void filter_nucleos(std::string& nucleos, const bool& upper) {
    if (upper) {
        for (char& c : nucleos) c = str_manip::upper_filter_table[c];
    } else {
        for (char& c : nucleos) c = str_manip::filter_table[c];
    }
    return;
}



// Trim leading and trailing whitespace off a string
inline void trimws(std::string& str) {
    size_t first = str.find_first_not_of(' ');
    if (std::string::npos == first) return;
    size_t last = str.find_last_not_of(' ');
    str = str.substr(first, (last - first + 1));
    return;
}


/*
 Split a string based on a single-character delimiter
 */

inline std::vector<std::string> cpp_str_split_delim(const std::string& in_string,
                                                    const char& split) {


    std::vector<std::string> out(1, "");
    std::string::size_type n = 1;

    std::string::size_type i = in_string.find(split);
    if (i != std::string::npos) {
        // Index for the output vector
        uint64 j = 0;
        // Index for the previous i:
        uint64 i0 = 0;
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

/*
 Split a string based on newlines. It also checks for Windows '\r\n' newlines and
 removes the '\r'.
 */
inline std::vector<std::string> cpp_str_split_newline(const std::string& in_string) {

    const char split = '\n';

    std::vector<std::string> out(1, "");
    std::string::size_type n = 1;

    std::string::size_type i = in_string.find(split);
    if (i != std::string::npos) {
        // Index for the output vector
        uint64 j = 0;
        // Index for the previous i:
        uint64 i0 = 0;
        while (i != std::string::npos) {
            for (std::string::size_type k = i0; k < i; k++) {
                out[j] += in_string[k];
            }
            if (out[j].back() == '\r') out[j].pop_back();
            i0 = i + n;
            i = in_string.find(split, i0);
            j++;
            out.push_back("");
        }
        for (std::string::size_type k = i0; k < in_string.size(); k++) {
            out[j] += in_string[k];
        }
        if (out[j].back() == '\r') out[j].pop_back();
    } else {
        out[0] = in_string;
        if (out[0].back() == '\r') out[0].pop_back();
    }

    return out;
}




// Split a main string by a string delimeter:
inline std::vector<std::string> cpp_str_split_delim_str(const std::string& in_string,
                                                        const std::string& split) {
    std::vector<std::string> splitted;
    size_t last = 0;
    size_t next = 0;
    while ((next = in_string.find(split, last)) != std::string::npos) {
        splitted.push_back(in_string.substr(last, next-last));
        last = next + split.size();
    }
    splitted.push_back(in_string.substr(last));
    return splitted;
}



// Count a substring in a main string:
inline uint32 count_substr(const std::string& in_string, const std::string& substr) {
    uint32 occurrences = 0;
    std::string::size_type pos = 0;
    while ((pos = in_string.find(substr, pos)) != std::string::npos) {
        ++ occurrences;
        pos += substr.length();
    }
    return occurrences;
}



/*
 Reverse complement of a DNA chromosome.

 Make sure that `chrom` contains only T, C, A, or G!
 */
inline void rev_comp(std::string& chrom) {

    uint64 n = chrom.size();
    uint64 half_n = n / 2;
    char tmp;

    for (uint64 j = 0; j < half_n; j++) {
        tmp = str_manip::cmp_map[chrom[j]]; // goes to `n-j-1`
        chrom[j] = str_manip::cmp_map[chrom[(n-j-1)]];
        chrom[(n-j-1)] = tmp;
    }

    if ((n & 1ULL) == 1ULL) chrom[half_n] = str_manip::cmp_map[chrom[half_n]];

    return;
}

/*
 Same thing, except that it only does it for the first `n` characters in `chrom`
 */
inline void rev_comp(std::string& chrom, const uint64& n) {

    uint64 half_n = n / 2;
    char tmp;

    for (uint64 j = 0; j < half_n; j++) {
        tmp = str_manip::cmp_map[chrom[j]]; // goes to `n-j-1`
        chrom[j] = str_manip::cmp_map[chrom[(n-j-1)]];
        chrom[(n-j-1)] = tmp;
    }

    if ((n & 1ULL) == 1ULL) chrom[half_n] = str_manip::cmp_map[chrom[half_n]];

    return;
}



# endif
