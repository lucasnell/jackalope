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

#include "gemino_types.h"


using namespace Rcpp;





// Needed for multiple functions below (it's faster than using std::pow(a, b))
// a^b
inline double fast_pow(double a, double b) {
    return std::exp(b * std::log(a));
}


namespace vitter {
    // Adding +1 to this so that when sampling, `engine() / seq_var::sitmo_max` is always
    // < 1. This is useful for discrete sampling bc I don't want there to be a number
    // that has a 1/2^32 chance of being sampled like what would happen if engine()
    // produced a number == `sitmo::prng_engine::max()`.
    const double sitmo_max = (double) sitmo::prng_engine::max() + 1.0;
}





// =============================================================================
// =============================================================================

//          Distribution of S

// =============================================================================
// =============================================================================



// Pr(S = s)
// (This is the same as `f` below, but vectorized and exported to R.
// I'm keeping them separate so I can keep `f` inline.)
//[[Rcpp::export]]
std::vector<double> f_s(const std::vector<double>& s_vec, const double& n, const double& N) {
    std::vector<double> out(s_vec.size());
    double out_i, s;
    for (uint i = 0; i < s_vec.size(); i++) {
        s = s_vec[i];
        if (s < 0 || s > (N - n)) {
            out[i] = 0;
            continue;
        }
        out_i = (n / N);
        for (int j = 0; j < s; j++) {
            out_i *= ((N - n - (double) j) / (N - 1 - (double) j));
        }
        out[i] = out_i;
    }
    return out;
}

// This is Pr(S <= s)
//[[Rcpp::export]]
std::vector<double> F_s(const std::vector<double>& s_vec, const double& n, const double& N) {
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





// var(S)
//[[Rcpp::export]]
double variance_s(double n, double N) {
    return ((N + 1) * (N - n) * n) / ((n + 2) * fast_pow(n + 1, 2));
}


// If comparing distances to S, remember that S in Vitter's paper is the number of
// positions to skip _before_ taking the next one, so it should be 1 less than the
// distance.
// So just add 1 to this function to get expected distances.

// Expected value of S
//[[Rcpp::export]]
double expected_s(double n, double N) {
    return (N - n) / (n + 1);
}


// Pr(S == s)
inline double f(const double& s, const double& n, const double& N) {
    if (s < 0 || s > (N - n)) return 0;
    double out = (n / N);
    for (int i = 0; i < s; i++) {
        out *= ((N - n - (double) i) / (N - 1 - (double) i));
    }
    return out;
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



// ============================
// Algorithm D_1
// (when n^2 / N <= the n2N parameter (default = 50))
// ============================

// --------
// Inline functions for Algorithm D_1
// --------
inline double g1(const double& x, const double& n, const double& N) {
    if (x < 0 || x > N) return 0;
    return (n / N) * fast_pow(1 - x/N, n - 1);
}
inline double c1(const double& n, const double& N) {
    return N / (N - n + 1);
}
inline double h1(const double& s, const double& n, const double& N) {
    if (s < 0 || s > (N - n)) return 0;
    return (n / N) * fast_pow(1 - (s / (N - n + 1)), n - 1);
}
inline double x1(const double& U, const double& n, const double& N) {
    return N * (1 - fast_pow(U, 1/n));
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





// ============================
// Algorithm D_2
// (when n^2 / N > the n2N parameter (default = 50))
// ============================


// --------
// Inline functions for Algorithm D_2
// --------
inline double g2(const double& s, const double& n, const double& N) {
    if (s < 0) stop("Computational error. s cannot < 0 in g2.");
    return ((n - 1) / (N - 1)) * fast_pow(1 - ((n - 1) / (N - 1)), s);
}
inline double c2(const double& n, const double& N) {
    return (n / (n - 1)) * ((N - 1) / N);
}
inline double h2(const double& s, const double& n, const double& N) {
    if (s < 0 || s > (N - n)) return 0;
    return (n / N) * fast_pow(1 - ((n - 1) / (N - s)), s);
}
inline double x2(const double& U, const double& n, const double& N) {
    return std::floor(std::log(U) / std::log(1 - ((n - 1) / (N - 1))));
}

// --------
// One S value for Algorithm D_2
// --------
uint algorithm_d2_S(const sint& n, const uint& N, sitmo::prng_engine& engine,
                    const double alpha) {

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






// ========================================================
// ========================================================

// Full Algorithm D

// ========================================================
// ========================================================


// From testing, alpha = 0.8 and n2N = 50 have fastest times
// this will return 0-->(N-1)

std::vector<uint> vitter_d(sint n, uint N, sitmo::prng_engine& engine,
                   double n2N = 50, double alpha = 0.8) {

    // Commented this out bc this will crash R if run in parallel and stop happens.
    // if (alpha > 1 || alpha < 0) stop("Invalid alpha. It must be (0,1).");

    std::vector<uint> samples;
    uint S;
    int_fast64_t current_pos = -1;
    if ((fast_pow(n, 2) / N) > n2N) {
        while (n > 1) {
            S = algorithm_d2_S(n, N, engine, alpha);
            current_pos += S + 1;
            samples.push_back(current_pos);
            n--;
            N -= (S + 1);
        }
        // At n = 1, D2 divides by zero, but this works just fine
        if (n == 1) {
            S = ((double) engine() / vitter::sitmo_max) * N;
            current_pos += S + 1;
            samples.push_back(current_pos);
        }
    } else {
        while (n > 0) {
            S = algorithm_d1_S(n, N, engine, alpha);
            current_pos += S + 1;
            samples.push_back(current_pos);
            n --;
            N -= (S + 1);
        }
    }

    return samples;
}


//[[Rcpp::export]]
arma::Mat<uint> test_vitter_d(const uint reps, uint n, uint N, const std::vector<uint> seeds,
                              const double n2N = 50, const double alpha = 0.8) {

    arma::Mat<uint> out(reps, n);
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
        std::vector<uint> point_positions = vitter_d(n, N, engine, n2N, alpha);
        if (point_positions.size() != n) stop("Inappropriate output.");
        for (uint j = 0; j < n; j++) {
            out(i,j) = point_positions[j];
        }
        // Rcpp::checkUserInterrupt(); // <-- Causes crash when in parallel
    }
    #ifdef _OPENMP
    }
    #endif

    return out;
}

