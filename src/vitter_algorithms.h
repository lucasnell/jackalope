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

template <class T>
void vitter_d(T& input_vec, uint N, sitmo::prng_engine& engine,
              const double n2N = 50, const double alpha = 0.8);


#endif
