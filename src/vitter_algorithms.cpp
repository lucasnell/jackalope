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
#include <pcg/pcg_random.hpp> // pcg prng

#ifdef _OPENMP
#include <omp.h>
#endif

#include "gemino_types.h"  // integer types
#include "util.h"  // uints_get_size, cpp_choose
#include "vitter_algorithms.h"  // vitter namespace
#include "pcg.h"


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
    double out = cpp_choose(N - s - 1.0, n - 1.0);
    out /= cpp_choose(N, n);
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
    for (uint32 i = 0; i < s_vec.size(); i++) {
        s = s_vec[i];
        if (s < 0 || s > (N - n)) {
            out[i] = 0;
            continue;
        }
        out_i = (N - n) / N;
        for (int j = 1; j <= s; j++) {
            out_i *= ((N - n - static_cast<double>(j)) / (N - static_cast<double>(j)));
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
uint32 vitter_a_S(double n, double N, pcg32& engine) {
    double V = runif_01(engine);
    uint32 s = 0;
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
    return (n / N) * std::pow(1 - x/N, static_cast<uint32>(n - 1));
}
inline double c1(const double& n, const double& N) {
    return N / (N - n + 1);
}
inline double h1(const double& s, const double& n, const double& N) {
    if (s < 0 || s > (N - n)) return 0;
    return (n / N) * std::pow(1 - (s / (N - n + 1)), static_cast<uint32>(n - 1));
}
inline double x1(const double& U, const double& n, const double& N) {
    return N * (1 - std::pow(U, 1/n));
}
// --------
// One S value for Algorithm D_1
// --------
uint32 algorithm_d1_S(const sint32& n, const uint32& N, pcg32& engine,
                    const double alpha) {

    double U, X, c, comp_denom;
    uint32 S;
    if (n < (alpha * N)) {
        while (true) {
            U = runif_01(engine);
            X = x1(runif_01(engine), static_cast<double>(n), static_cast<double>(N));
            c = c1(static_cast<double>(n), static_cast<double>(N));
            comp_denom = c * g1(X, static_cast<double>(n), static_cast<double>(N));
            if (U <= h1(std::floor(X), static_cast<double>(n),
                        static_cast<double>(N)) / comp_denom) {
                S = std::floor(X);
                break;
            } else if (U <= f(std::floor(X), static_cast<double>(n),
                              static_cast<double>(N)) / comp_denom) {
                S = std::floor(X);
                break;
            } else {
                continue;
            }
        }
    } else {
        S = vitter_a_S(static_cast<double>(n), static_cast<double>(N), engine);
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
    return ((n - 1) / (N - 1)) * std::pow(1 - ((n - 1) / (N - 1)), static_cast<uint32>(s));
}
inline double c2(const double& n, const double& N) {
    return (n / (n - 1)) * ((N - 1) / N);
}
inline double h2(const double& s, const double& n, const double& N) {
    if (s < 0 || s > (N - n)) return 0;
    return (n / N) * std::pow(1 - ((n - 1) / (N - s)), static_cast<uint32>(s));
}
inline double x2(const double& U, const double& n, const double& N) {
    return std::floor(std::log(U) / std::log(1 - ((n - 1) / (N - 1))));
}

// --------
// One S value for Algorithm D_2
// --------
uint32 algorithm_d2_S(const sint32& n, const uint32& N, pcg32& engine,
                    const double& alpha) {

    double U, X, c, comp_denom;
    uint32 S;
    if (n < (alpha * N)) {
        while (true) {
            U = runif_01(engine);
            X = x2(runif_01(engine), static_cast<double>(n), static_cast<double>(N));
            c = c2(static_cast<double>(n), static_cast<double>(N));
            comp_denom = c * g2(X, static_cast<double>(n), static_cast<double>(N));
            if (U <= h2(std::floor(X), static_cast<double>(n),
                        static_cast<double>(N)) / comp_denom) {
                S = std::floor(X);
                break;
            } else if (U <= f(std::floor(X), static_cast<double>(n),
                              static_cast<double>(N)) / comp_denom) {
                S = std::floor(X);
                break;
            } else {
                continue;
            }
        }
    } else {
        S = vitter_a_S(static_cast<double>(n), static_cast<double>(N), engine);
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
arma::Mat<uint32> test_vitter_d(const uint32 reps, uint32 n, uint32 N,
                                const uint32& n_cores,
                                const double& n2N = 50, const double& alpha = 0.8) {

    arma::Mat<uint32> out(n, reps);
    if (alpha > 1 || alpha < 0) stop("Invalid alpha. It must be (0,1).");
    if (n > N) stop("n must be <= N.");

    const std::vector<std::vector<uint64>> seeds = mc_seeds(n_cores);

    #ifdef _OPENMP
    #pragma omp parallel num_threads(n_cores) if (n_cores > 1)
    {
    #endif

    std::vector<uint64> active_seeds;

    // Write the active seed per core or just write one of the seeds.
    #ifdef _OPENMP
    active_seeds = seeds[omp_get_thread_num()];
    #else
    active_seeds = seeds[0];
    #endif

    pcg32 engine = seeded_pcg(active_seeds);

    // Parallelize the Loop
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint32 i = 0; i < reps; i++) {
        arma::uvec point_positions(n);
        vitter_d<arma::uvec>(point_positions, N, engine, n2N, alpha);
        for (uint32 j = 0; j < n; j++) out(j,i) = static_cast<uint32>(point_positions(j));
        // Rcpp::checkUserInterrupt(); // <-- Causes crash when in parallel
    }
    #ifdef _OPENMP
    }
    #endif

    return out;
}

