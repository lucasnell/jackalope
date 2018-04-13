# ifndef __GEMINO_STR_MANIP_H
# define __GEMINO_STR_MANIP_H


#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

using namespace Rcpp;


void cpp_to_upper(std::string& input_str);

std::string cpp_merge_str(const std::vector<std::string>& in_strings);

std::vector<std::string> cpp_str_split_delim(const std::string& in_string,
                                             const char& split);


# endif
