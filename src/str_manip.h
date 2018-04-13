# ifndef __GEMINO_STR_MANIP_H
# define __GEMINO_STR_MANIP_H


#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

using namespace Rcpp;


std::string to_upper(const std::string& input_str);
void cpp_to_upper(std::string& input_str);


std::string uvec_to_str(const std::vector<uint8_t>& in_uvec);
std::string ucol_to_str(const arma::Col<uint8_t>& in_ucol);

std::vector<uint8_t> str_to_uvec(const std::string& in_str);
arma::Col<uint8_t> str_to_ucol(const std::string& in_str);


std::string cpp_merge_str(const std::vector<std::string>& in_strings);

std::vector<std::string> cpp_str_split_int(const std::string& in_string, const int& n = 1);


std::vector<std::string> cpp_str_split_delim(const std::string& in_string, const char& split);

# endif
