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


/*
 ============================
 Algorithm D_1
 (when n^2 / N <= the n2N parameter (default = 50))
 ============================
 */
uint algorithm_d1_S(const sint& n, const uint& N, sitmo::prng_engine& engine,
                    const double alpha);

/*
 ============================
 Algorithm D_2
 (when n^2 / N > the n2N parameter (default = 50))
 ============================
 */
uint algorithm_d2_S(const sint& n, const uint& N, sitmo::prng_engine& engine,
                    const double& alpha);

/*
 ============================
 Full algorithm
 ============================
*/
template <class T>
void vitter_d(T& input_vec, uint N, sitmo::prng_engine& engine,
              const double n2N = 50, const double alpha = 0.8);


#endif
