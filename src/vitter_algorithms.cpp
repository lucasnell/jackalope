//
// Algorithms for fast sampling without replacement
//
// These algorithms are from the following paper:
// Vitter, Jeffrey Scott. 1984. Faster methods for random sampling. Communications of
//     the ACM 27:703–718.
//
// My goal was to get Algorithm D working, so I haven't included full versions of
// any others.
//
// Here are the variables referred to in the code below:
//
// S: The number of positions to skip before sampling the next available position.
//     For example, if S = 0 and the current position is 100, the next position
//     sampled will be 101.
// N: The population size. The sampling will generate numbers from 0 to (N - 1).
// n: The sample size. Sampling will generate n random numbers.
//     n should always be <= N.


#include <RcppArmadillo.h>
#include <vector>
#include <numeric>
#include <cmath>
#include <sitmo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "gemino_types.h"  // integer types
#include "util.h"  // uints_get_size
#include "vitter_algorithms.h"  // vitter namespace


using namespace Rcpp;








// =============================================================================
// =============================================================================

//          Distribution of S

// =============================================================================
// =============================================================================



//' Pr(S == s).
//'
//'
//' @noRd
inline double f(const double& s, const double& n, const double& N) {
    if (s < 0 || s > (N - n)) return 0;
    double out = (n / N);
    for (int i = 0; i < s; i++) {
        out *= ((N - n - (double) i) / (N - 1 - (double) i));
    }
    return out;
}


//' Pr(S <= s)
//'
//' @noRd
//'
std::vector<double> F_s(const std::vector<double>& s_vec, const double& n,
                        const double& N) {
    std::vector<double> out(s_vec.size());
    double out_i, s;
    for (uint i = 0; i < s_vec.size(); i++) {
        s = s_vec[i];
        if (s < 0 || s > (N - n)) {
            out[i] = 0;
            continue;
        }
        out_i = (N - n) / N;
        for (int j = 1; j <= s; j++) {
            out_i *= ((N - n - (double) j) / (N - (double) j));
        }
        out[i] = 1 - out_i;
    }
    return out;
}





//' var(S).
//'
//' @noRd
//'
double variance_s(double n, double N) {
    return ((N + 1) * (N - n) * n) / ((n + 2) * ((n + 1) * (n + 1)));
}

//' Expected value of S.
//'
//' If comparing distances to S, remember that S in Vitter's paper is the number of
//' positions to skip _before_ taking the next one, so it should be 1 less than the
//' distance.
//' So just add 1 to this function to get expected distances.
//'
//' @noRd
//'
double expected_s(double n, double N) {
    return (N - n) / (n + 1);
}








// --------
// One S value using Algorithm A (used in Algorithm D if n >= alpha * N)
// --------
uint vitter_a_S(double n, double N, sitmo::prng_engine& engine) {
    double V = (double) engine() / vitter::sitmo_max;
    uint s = 0;
    double lhs = N - n;
    double rhs = V * N;
    while (lhs > rhs) {
        s++;
        lhs *= (N - n - s);
        rhs *= (N - s);
    }
    return s;
}



/*
 ============================
 Algorithm D_1
 (when n^2 / N <= the n2N parameter (default = 50))
 ============================
 */

// --------
// Inline functions for Algorithm D_1
// --------
inline double g1(const double& x, const double& n, const double& N) {
    if (x < 0 || x > N) return 0;
    return (n / N) * std::pow(1 - x/N, n - 1);
}
inline double c1(const double& n, const double& N) {
    return N / (N - n + 1);
}
inline double h1(const double& s, const double& n, const double& N) {
    if (s < 0 || s > (N - n)) return 0;
    return (n / N) * std::pow(1 - (s / (N - n + 1)), n - 1);
}
inline double x1(const double& U, const double& n, const double& N) {
    return N * (1 - std::pow(U, 1/n));
}
// --------
// One S value for Algorithm D_1
// --------
uint algorithm_d1_S(const sint& n, const uint& N, sitmo::prng_engine& engine,
                    const double alpha) {

    double U, X, c, comp_denom;
    uint S;
    if (n < (alpha * N)) {
        while (true) {
            U = (double) engine() / vitter::sitmo_max;
            X = x1((double) engine() / vitter::sitmo_max, (double) n, (double) N);
            c = c1((double) n, (double) N);
            comp_denom = c * g1(X, (double) n, (double) N);
            if (U <= h1(std::floor(X), (double) n, (double) N) / comp_denom) {
                S = std::floor(X);
                break;
            } else if (U <= f(std::floor(X), (double) n, (double) N) / comp_denom) {
                S = std::floor(X);
                break;
            } else {
                continue;
            }
        }
    } else {
        S = vitter_a_S((double) n, (double) N, engine);
    }
    return S;
}





/*
 ============================
 Algorithm D_2
 (when n^2 / N > the n2N parameter (default = 50))
  ============================
 */


// --------
// Inline functions for Algorithm D_2
// --------
inline double g2(const double& s, const double& n, const double& N) {
    if (s < 0) stop("Computational error. s cannot < 0 in g2.");
    return ((n - 1) / (N - 1)) * std::pow(1 - ((n - 1) / (N - 1)), s);
}
inline double c2(const double& n, const double& N) {
    return (n / (n - 1)) * ((N - 1) / N);
}
inline double h2(const double& s, const double& n, const double& N) {
    if (s < 0 || s > (N - n)) return 0;
    return (n / N) * std::pow(1 - ((n - 1) / (N - s)), s);
}
inline double x2(const double& U, const double& n, const double& N) {
    return std::floor(std::log(U) / std::log(1 - ((n - 1) / (N - 1))));
}

// --------
// One S value for Algorithm D_2
// --------
uint algorithm_d2_S(const sint& n, const uint& N, sitmo::prng_engine& engine,
                    const double& alpha) {

    double U, X, c, comp_denom;
    uint S;
    if (n < (alpha * N)) {
        while (true) {
            U = (double) engine() / vitter::sitmo_max;
            X = x2((double) engine() / vitter::sitmo_max, (double) n, (double) N);
            c = c2((double) n, (double) N);
            comp_denom = c * g2(X, (double) n, (double) N);
            if (U <= h2(std::floor(X), (double) n, (double) N) / comp_denom) {
                S = std::floor(X);
                break;
            } else if (U <= f(std::floor(X), (double) n, (double) N) / comp_denom) {
                S = std::floor(X);
                break;
            } else {
                continue;
            }
        }
    } else {
        S = vitter_a_S((double) n, (double) N, engine);
    }
    return S;
}





/*
 ========================================================
 ========================================================

 Full Algorithm D

 ========================================================
 ========================================================
 */





//[[Rcpp::export]]
arma::Mat<uint> test_vitter_d(const uint reps, uint n, uint N,
                              const std::vector<uint> seeds,
                              const double n2N = 50, const double alpha = 0.8) {

    arma::Mat<uint> out(n, reps);
    if (alpha > 1 || alpha < 0) stop("Invalid alpha. It must be (0,1).");
    if (n > N) stop("n must be <= N.");

    const uint n_cores = seeds.size();

    #ifdef _OPENMP
    #pragma omp parallel num_threads(n_cores) if (n_cores > 1)
    {
    #endif

    uint active_seed;

    // Write the active seed per core or just write one of the seeds.
    #ifdef _OPENMP
    active_seed = seeds[omp_get_thread_num()];
    #else
    active_seed = seeds[0];
    #endif

    sitmo::prng_engine engine(active_seed);

    // Parallelize the Loop
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint i = 0; i < reps; i++) {
        arma::uvec point_positions(n);
        vitter_d<arma::uvec>(point_positions, N, engine, n2N, alpha);
        out.col(i) = point_positions;
        // Rcpp::checkUserInterrupt(); // <-- Causes crash when in parallel
    }
    #ifdef _OPENMP
    }
    #endif

    return out;
}

