
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>



using namespace Rcpp;


template<typename T_it>
void SequenceToUpperCase( T_it begin, T_it end )
{
    // Convert to upper: clear the '32' bit, 0x20 in hex. And with the
    // inverted bit string (~).
    for ( auto it = begin; it != end; ++it )
        *it &= ~0x20;
}
//' Make a string uppercase in place.
//'
//' This is ~4x faster than R's .Internal(toupper(...))
//'
//' @param input_str An input string.
//'
//' @return Output string of all uppercase letters.
//'
//[[Rcpp::export]]
std::string to_upper(const std::string& input_str) {
    std::string output_str = input_str;
    SequenceToUpperCase(output_str.begin(), output_str.end());
    return output_str;
}

// Same thing, but for C++ only
void cpp_to_upper(std::string& input_str) {
    SequenceToUpperCase(input_str.begin(), input_str.end());
    return;
}



//' Make a string from a vector of integers.
//'
//'
//'
//' @param in_vec A vector of unsigned, 8-bit integers.
//'
//' @return A string.
//'
//'
// [[Rcpp::export]]
std::string uvec_to_str(const std::vector<uint8_t>& in_uvec) {
    std::string out_str = "";
    for (int i = 0; i < in_uvec.size(); i++) {
        out_str += in_uvec[i];
    }
    return out_str;
}


//' Make a string from a 1-column matrix of integers.
//'
//'
//' @param in_vec A 1-column matrix of unsigned, 8-bit integers.
//'
//' @return A string.
//'
//'
// [[Rcpp::export]]
std::string ucol_to_str(const arma::Col<uint8_t>& in_ucol) {
    std::string out_str = "";
    for (int i = 0; i < in_ucol.size(); i++) {
        out_str += in_ucol(i);
    }
    return out_str;
}



//' Make a vector of integers from a string.
//'
//'
//' @param in_str A string.
//'
//' @return A vector of unsigned, 8-bit integers.
//'
//'
// [[Rcpp::export]]
std::vector<uint8_t> str_to_uvec(const std::string& in_str) {
    std::vector<uint8_t> out_uvec(in_str.size());
    for (int i = 0; i < in_str.size(); i++) {
        out_uvec[i] = in_str[i];
    }
    return out_uvec;
}

//' Make a 1-column matrix of integers from a string.
//'
//'
//' @param in_str A string.
//'
//' @return A 1-column matrix of unsigned, 8-bit integers.
//'
//'
// [[Rcpp::export]]
arma::Col<uint8_t> str_to_ucol(const std::string& in_str) {
    arma::Col<uint8_t> out_ucol(in_str.size());
    for (int i = 0; i < in_str.size(); i++) {
        out_ucol[i] = in_str[i];
    }
    return out_ucol;
}





//' Merge a vector of strings into one.
//'
//'
//' @param in_strings Character vector of strings to merge.
//'
//' @return A single string.
//'
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






//' Split one string into a vector of strings of the same size.
//'
//' \emph{Note}: This function is not exported to R.
//'
//' \emph{Note}: The last string may be of a different length if
//'     \code{in_string.length() / n != round(in_string.length() / n)}.
//'
//' @param in_string Input string to split.
//' @param n Size of strings to output.
//'
//' @return Vector of strings of the same size.
//'
//' @name cpp_str_split_int
//'
std::vector<std::string> cpp_str_split_int(const std::string& in_string, const int& n = 1) {

    int num_substr = in_string.length() / n;

    std::vector<std::string> out(num_substr);

    for (int j = 0; j < num_substr; j++) {
        out[j] = in_string.substr(j * n, n);
    }

    return out;
}



// Split a string based on a single-character delimiter

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
