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

#include "gemino_types.hpp"
#include "util.hpp"  // uints_get_size
#include "pcg.hpp"  // pcg::max


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






// --------
// One S value using Algorithm A (used in Algorithm D if n >= alpha * N)
// --------
inline uint32 vitter_a_S(double n, double N, pcg64& engine) {
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
inline uint32 algorithm_d1_S(const sint32& n, const uint32& N, pcg64& engine,
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
inline uint32 algorithm_d2_S(const sint32& n, const uint32& N, pcg64& engine,
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
//'     `std::vector<uint32>`) of length `n`.
//'     Sampling will generate `n` random numbers. `n` should always be <= N.
//' @param N The population size. The sampling will generate numbers from
//'     `0` to `(N - 1)`.
//' @param engine A pcg PRNG engine.
//' @param start An unsigned integer to add to each position. In effect, this makes
//'     the output a range from `start` to `start + N - 1`.
//'     Defaults to `0`.
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
void vitter_d(T& samples, uint32 N, pcg64& engine,
              const uint32& start = 0,
              const double n2N = 50, const double alpha = 0.8) {

    // Commented this out bc this will crash R if run in parallel and stop happens.
    // if (alpha > 1 || alpha < 0) stop("Invalid alpha. It must be (0,1).");

    sint32 n = static_cast<sint32>(uints_get_size(samples));

    uint32 S, ind = 0;
    sint64 current_pos = -1;
    if (((n * n) / N) > n2N) {
        while (n > 1) {
            S = algorithm_d2_S(n, N, engine, alpha);
            current_pos += S + 1;
            samples[ind] = current_pos + start;
            ind++;
            n--;
            N -= (S + 1);
        }
        // At n = 1, D2 divides by zero, but this works just fine
        if (n == 1) {
            S = runif_01(engine) * N;
            current_pos += S + 1;
            samples[ind] = current_pos + start;
            ind++;
        }
    } else {
        while (n > 0) {
            S = algorithm_d1_S(n, N, engine, alpha);
            current_pos += S + 1;
            samples[ind] = current_pos + start;
            ind++;
            n--;
            N -= (S + 1);
        }
    }

    return;
}



#endif
