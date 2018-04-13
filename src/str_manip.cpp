
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>



using namespace Rcpp;




//' Make a string uppercase in place.
//'
//' The version of this function exported to R was ~4x faster than R's
//' `.Internal(toupper(...))`, and keeping it in C++ is surely even faster.
//'
//' @param input_str An input string to be changed.
//'
//' @return Nothing. Changes are made in place.
//'
//' @noRd
//'
void cpp_to_upper(std::string& input_str) {
    /*
     Convert to upper: clear the '32' bit, 0x20 in hex. And with the
     inverted bit string (~).
     */
    for (std::string::iterator it = input_str.begin(); it != input_str.begin(); ++it) {
        *it &= ~0x20;
    }
    return;
}





//' Merge a vector of strings into one.
//'
//'
//' @param in_strings Character vector of strings to merge.
//'
//' @return A single string.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::string cpp_merge_str(const std::vector<std::string>& in_strings) {

    int num_char = in_strings.size();
    std::string out_str;

    for(int j=0; j < num_char; j++) {
        out_str += in_strings[j];
    }

    return out_str;
}





//' Split a string based on a single-character delimiter
//'
//'
//' @param in_string A string to split.
//' @param split Character to split string by.
//'
//' @noRd
//'
//[[Rcpp::export]]
std::vector<std::string> cpp_str_split_delim(const std::string& in_string,
                                             const char& split) {


    std::vector<std::string> out(1, "");
    std::string::size_type n = 1;

    std::string::size_type i = in_string.find(split);
    if (i != std::string::npos) {
        // Index for the output vector
        uint j = 0;
        // Index for the previous i:
        uint i0 = 0;
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
