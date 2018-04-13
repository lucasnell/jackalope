# ifndef __GEMINO_VITTER_ALGS_H
# define __GEMINO_VITTER_ALGS_H

#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <numeric>
#include <cmath>
#include <sitmo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "gemino_types.h"


using namespace Rcpp;




std::vector<double> f_s(const std::vector<double>& s_vec, const double& n, const double& N);
std::vector<double> F_s(const std::vector<double>& s_vec, const double& n, const double& N);
double variance_s(double n, double N);
double expected_s(double n, double N);

uint vitter_a_S(double n, double N, sitmo::prng_engine& engine);
uint algorithm_d1_S(const sint& n, const uint& N, sitmo::prng_engine& engine,
                    const double alpha);
uint algorithm_d2_S(const sint& n, const uint& N, sitmo::prng_engine& engine,
                    const double alpha);
std::vector<uint> vitter_d(sint n, uint N, sitmo::prng_engine& engine,
                           double n2N = 50, double alpha = 0.8);


#endif
