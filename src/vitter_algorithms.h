# ifndef __GEMINO_VITTER_ALGS_H
# define __GEMINO_VITTER_ALGS_H

/*
 Algorithms for fast sampling *without* replacement
 */

#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <numeric>
#include <cmath>
#include <pcg/pcg_random.hpp> // pcg prng

#ifdef _OPENMP
#include <omp.h>
#endif

#include "gemino_types.h"
#include "util.h"  // uints_get_size
#include "pcg.h"  // pcg::max


using namespace Rcpp;




/*
 ============================
 Algorithm D_1
 (when n^2 / N <= the n2N parameter (default = 50))
 ============================
 */
uint algorithm_d1_S(const sint& n, const uint& N, pcg32& engine,
                    const double alpha);

/*
 ============================
 Algorithm D_2
 (when n^2 / N > the n2N parameter (default = 50))
 ============================
 */
uint algorithm_d2_S(const sint& n, const uint& N, pcg32& engine,
                    const double& alpha);

/*
 ============================
 Full algorithm
 ============================
*/
//' "Algorithm D" for fast sampling without replacement.
//'
//' This algorithm is from the following paper:
//' Vitter, Jeffrey Scott. 1984. Faster methods for random sampling. Communications of
//'     the ACM 27:703–718.
//'
//' @param input_vec A vector of unsigned integers (class `arma::uvec` or
//'     `std::vector<uint>`) of length `n`.
//'     Sampling will generate `n` random numbers. `n` should always be <= N.
//' @param N The population size. The sampling will generate numbers from
//'     `0` to `(N - 1)`.
//' @param engine A pcg PRNG engine.
//' @param n2N A numeric threshold placed on the algorithm used to find new locations.
//'     This is not recommended to be changed. Defaults to 50.
//' @param alpha A numeric threshold placed on the algorithm used to find new locations.
//'     This is not recommended to be changed. Defaults to 0.8.
//'
//'
//'
//'
//' @noRd
//'
//'
template <typename T>
void vitter_d(T& samples, uint N, pcg32& engine,
              const double n2N = 50, const double alpha = 0.8) {

    // Commented this out bc this will crash R if run in parallel and stop happens.
    // if (alpha > 1 || alpha < 0) stop("Invalid alpha. It must be (0,1).");

    sint n = static_cast<sint>(uints_get_size(samples));

    uint S, ind = 0;
    sint64 current_pos = -1;
    if (((n * n) / N) > n2N) {
        while (n > 1) {
            S = algorithm_d2_S(n, N, engine, alpha);
            current_pos += S + 1;
            samples[ind] = current_pos;
            ind++;
            n--;
            N -= (S + 1);
        }
        // At n = 1, D2 divides by zero, but this works just fine
        if (n == 1) {
            S = runif_01(engine) * N;
            current_pos += S + 1;
            samples[ind] = current_pos;
            ind++;
        }
    } else {
        while (n > 0) {
            S = algorithm_d1_S(n, N, engine, alpha);
            current_pos += S + 1;
            samples[ind] = current_pos;
            ind++;
            n--;
            N -= (S + 1);
        }
    }

    return;
}



#endif
