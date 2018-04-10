# ifndef _UTIL_H
# define _UTIL_H

#include <RcppArmadillo.h>
#include <vector>
#include <string>

#include "gemino_types.h"
// #include "new_variants.h" // new classes
// #include "alias.h" // alias sampling


using namespace Rcpp;


arma::ivec sample_ivec(arma::ivec in_vector, int size, bool replace,
                          arma::vec prob = arma::zeros<arma::vec>(0));

arma::uvec sample_uvec(arma::uvec in_vector, int size, bool replace,
                       arma::vec prob = arma::zeros<arma::vec>(0));

std::vector<uint> sample_big_ints(
        std::vector<uint> in_vector,
        bool replace = false,
        std::vector<uint> prob = std::vector<uint>(0));

int sample_int(arma::ivec in_vector, bool replace = false,
               arma::vec prob = arma::zeros<arma::vec>(0));


std::vector<int> str_to_intvec(const std::string& s);


std::string cpp_rando_seq(uint len);


double gc_prop(const std::string& sequence);
double gc_prop(const std::string& sequence, const uint& start, const uint& stop);

# endif
