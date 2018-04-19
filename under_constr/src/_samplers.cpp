#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <algorithm>  // lower_bound

#include "pcg/pcg_extras.hpp" // pcg prng
#include "pcg/pcg_random.hpp" // pcg prng


//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(sitmo, RcppArmadillo, gemino)]]

using namespace Rcpp;


#define SMALL_TOLERANCE 0.000000000001
#define MAX_UINT 4294967295

typedef uint_fast16_t uint16;
typedef uint_fast32_t uint;
typedef uint_fast64_t uint64;
typedef pcg_extras::pcg128_t uint128;

namespace samplers {
    double pcg_max = static_cast<double>(pcg32::max());
    long double pcg_max64 = static_cast<double>(pcg64::max());
    long double max64 = 18446744073709551616.0L;
}


/*
 THIS VERSION HAS WORKING 64-BIT SAMPLERS, BUT THE `pcg64_fast` DOESN'T WORK FOR SOME
 REASON. TESTING USING THE FIRST TWO FUNCTIONS HERE SHOW THAT THIS RNG DOES WORK, BUT
 I'M NOT IMPLEMENTING IT CORRECTLY IN SOME WAY.
 THE `pcg64_fast` DOESN'T APPEAR TO BE FASTER ANYWAY, SO I'M ABANDONDING IT.
 */



//[[Rcpp::export]]
std::vector<long double> test_fast64(const uint64& n, const std::vector<long double>& seeds) {

    std::vector<double> out(n);

    uint128 seed1 = samplers::pcg_max64 * seeds[0];
    uint128 seed2 = samplers::pcg_max64 * seeds[1];
    uint128 seed = (seed1<<64) + seed2;

    pcg64_fast eng(seed);

    std::vector<long double> rnd(n);

    for (uint64 i = 0; i < n; i++) {
        rnd[i] = static_cast<long double>(eng()) / samplers::pcg_max64;
    }

    return rnd;
}


//[[Rcpp::export]]
std::vector<long double> test_pcg64(const uint64& n, const std::vector<long double>& seeds) {

    uint128 seed1 = samplers::pcg_max64 * seeds[0];
    uint128 seed2 = samplers::pcg_max64 * seeds[1];
    uint128 seed = (seed1<<64) + seed2;

    pcg64 eng(seed);

    std::vector<long double> rnd(n);

    for (uint64 i = 0; i < n; i++) {
        rnd[i] = static_cast<long double>(eng()) / samplers::pcg_max64;
    }

    return rnd;
}


// Uses 4 calls to R::unif_rand to make two 64-bit seeds for a pcg32 RNG
pcg32 seeded_pcg() {

    // 32-bit seeds from unif_rand
    std::vector<uint64> sub_seeds = as<std::vector<uint64>>(Rcpp::runif(4,0,4294967296));
    uint64 seed1 = (sub_seeds[0]<<32) + sub_seeds[1];
    uint64 seed2 = (sub_seeds[2]<<32) + sub_seeds[3];

    pcg32 out(seed1, seed2);
    return out;
}

template <typename T>
uint sample_rare_(SEXP xptr_sexp, const uint64& N, const uint& rare) {

    XPtr<T> xptr(xptr_sexp);

    uint rares = 0;

    pcg32 eng = seeded_pcg();

    for (uint64 i = 0; i < N; i++) {
        uint k = xptr->sample(eng);
        if (k == rare) rares++;
    }

    return rares;
}



/*
 -----------
 64-bit versions
 -----------
 */

// Uses 8 calls to R::unif_rand to make two 128-bit seeds for a pcg64 RNG
template <typename RNG>
RNG seeded_pcg64() {

    // 32-bit seeds from unif_rand
    std::vector<uint64> sub_seeds = as<std::vector<uint64>>(Rcpp::runif(8,0,4294967296));
    uint64 seed64_1 = (sub_seeds[0]<<32) + sub_seeds[1];
    uint64 seed64_2 = (sub_seeds[2]<<32) + sub_seeds[3];
    // Commented these out so I can try this template with pcg64_fast, which only uses
    // one seed
    // uint64 seed64_3 = (sub_seeds[4]<<32) + sub_seeds[5];
    // uint64 seed64_4 = (sub_seeds[6]<<32) + sub_seeds[7];

    uint128 seed1 = (static_cast<uint128>(seed64_1)<<64) + static_cast<uint128>(seed64_2);
    // uint128 seed2 = (static_cast<uint128>(seed64_3)<<64) + static_cast<uint128>(seed64_4);

    RNG out(seed1);
    return out;
}










/*
 ========================================================================================
 ========================================================================================

 Table sampling from...
 Marsaglia, G., W. W. Tsang, and J. Wang. 2004. Fast generation of discrete random
 variables. Journal of Statistical Software 11.

 ========================================================================================
 ========================================================================================
 */


/*
 Converts a `p` vector of probabilities to a `ints` vector of integers, where each
 integer represents the approximate expected value of "successes" from 2^32 runs.
 If the sum of `ints` is != 2^32, this randomly chooses values to change based on
 probabilities in `p`.
 I did it this way because larger probabilities should be less affected by changing
 their respective values in `ints`.
 */
void fill_ints(const std::vector<double>& p, std::vector<uint>& ints) {

    // Vector holding the transitory values that will eventually be inserted into `int`
    arma::vec pp(p);
    pp /= arma::accu(pp);
    pp *= static_cast<double>(1L<<32);
    pp = arma::round(pp);

    // Converting to `ints`
    ints = arma::conv_to<std::vector<uint>>::from(pp);

    double d = static_cast<double>(1L<<32) - arma::accu(pp);

    // Vector for weighted sampling from vector of probabilities
    arma::vec p2(p);
    p2 /= arma::accu(p2);
    /*
     I'm not going to sample rare probabilities so anything < 2^-8 is set to zero
     for this sampling
     */
    double z = 1 / std::pow(2, 8);
    arma::uvec iv = arma::find(p2 < z);
    /*
     If there aren't any *above* this threshold, keep adding 8 to `x` in the expression
     `2^-x` until we would no longer be setting all probabilities to zero
     */
    while (iv.n_elem == p2.n_elem) {
        for (uint zz = 0; zz < 8; zz++) z *= 0.5;
        iv = arma::find(p2 < z);
    }
    p2(iv).fill(0);
    p2 /= arma::accu(p2);
    p2 = arma::cumsum(p2);

    // We need to remove from `ints`
    while (d < 0) {
        double u = R::unif_rand();
        iv = arma::find(p2 >= u, 1);
        ints[iv(0)]--;
        d++;
    }
    // We need to add to `ints`
    while (d > 0) {
        double u = R::unif_rand();
        iv = arma::find(p2 >= u, 1);
        ints[iv(0)]++;
        d--;
    }

    return;
}

class TableTable {
public:
    // Stores vectors of each category's Pr(sampled):
    std::vector<std::vector<uint>> T;
    // Stores values at which to transitiion between vectors of `T`:
    std::vector<uint> t;

    TableTable(const std::vector<double>& probs) : T(4), t(3, 0) {

        uint n_tables = T.size();

        uint n = probs.size();
        std::vector<uint> ints(n);
        // Filling the `ints` vector based on `probs`
        fill_ints(probs, ints);

        std::vector<uint> sizes(n_tables, 0);
        // Adding up sizes of `T` vectors:
        for (uint i = 0; i < n; i++) {
            for (uint k = 1; k <= n_tables; k++) {
                sizes[k-1] += dg(ints[i], k);
            }
        }
        // Adding up thresholds in the `t` vector
        for (uint k = 0; k < (n_tables - 1); k++) {
            t[k] = sizes[k]<<(32-8*(1+k));
            if (k > 0) t[k] += t[k-1];
        }
        // Re-sizing `T` vectors:
        for (uint i = 0; i < n_tables; i++) T[i].resize(sizes[i]);

        // Filling `T` vectors
        for (uint k = 1; k <= n_tables; k++) {
            uint ind = 0; // index inside `T[k-1]`
            for (uint i = 0; i < n; i++) {
                uint z = dg(ints[i], k);
                for (uint j = 0; j < z; j++) T[k-1][ind + j] = i;
                ind += z;
            }
        }
    }

    uint sample(pcg32& eng) {
        uint j = eng();
        if (j<t[0]) return T[0][j>>24];
        if (j<t[1]) return T[1][(j-t[0])>>(32-8*2)];
        if (j<t[2]) return T[2][(j-t[1])>>(32-8*3)];
        // for (uint i = 1; i < (T.size() - 1); i++) {
        //     if (j<t[i]) return T[i][(j-t[i-1])>>(32-8*(i+1))];
        // }
        return T[3][j-t[2]];
    }

    void print() const {
        // names coincide with names from Marsaglia (2004)
        std::vector<std::string> names = {"AA", "BB", "CC", "DD"};
        for (uint i = 0; i < T.size(); i++) {
            arma::urowvec x(T[i]);
            x.print(names[i] + ":");
        }
        arma::urowvec x(t);
        x.print("t:");
    }

private:
    static uint dg(const uint& m, const uint& k) {
        uint x = ((m>>(32-8*k))&255);
        return x;
    }
};



/*
 ---------------------------------------------------------------------------
 ---------------------------------------------------------------------------
 64-bit version
 ---------------------------------------------------------------------------
 ---------------------------------------------------------------------------
 */
// Check for less than in a long double vector
inline std::vector<uint> ld_lessthan(const std::vector<long double>& v,
                                     const long double& x) {
    std::vector<uint> out;
    for (uint i = 0; i < v.size(); i++) {
        if (v[i] < x) out.push_back(i);
    }
    return out;
}

template <typename RNG>
void fill_ints64(const std::vector<long double>& p, std::vector<uint64>& ints,
                 RNG& eng) {

    // Vector holding the transitory values that will eventually be inserted into `int`
    std::vector<long double> pp(p);
    ints.resize(pp.size());
    long double pp_sum = std::accumulate(pp.begin(), pp.end(), 0.0L);
    for (uint i = 0; i < pp.size(); i++) {
        pp[i] /= pp_sum;
        pp[i] *= samplers::max64;
        pp[i] = std::round(pp[i]);
        ints[i] = static_cast<uint64>(pp[i]);
    }
    pp_sum = std::accumulate(pp.begin(), pp.end(), 0.0L);

    long double d = samplers::max64 - pp_sum;

    // Vector for weighted sampling from vector of probabilities
    std::vector<long double> p2(p);
    long double p2_sum = std::accumulate(p2.begin(), p2.end(), 0.0L);
    for (uint i = 0; i < p2.size(); i++) p2[i] /= p2_sum;
    /*
    I'm not going to sample rare probabilities so anything < 2^-8 is set to zero
    for this sampling
    */
    long double z = 1 / std::pow(2, 8);
    std::vector<uint> iv = ld_lessthan(p2, z);
    /*
    If there aren't any *above* this threshold, keep adding 8 to `x` in the expression
    `2^-x` until we would no longer be setting all probabilities to zero
    */
    while (iv.size() == p2.size()) {
        for (uint zz = 0; zz < 8; zz++) z *= 0.5L;
        iv = ld_lessthan(p2, z);
    }
    for (uint i : iv) p2[i] = 0.0L;
    p2_sum = std::accumulate(p2.begin(), p2.end(), 0.0L);
    for (uint i = 0; i < p2.size(); i++) p2[i] /= p2_sum;
    for (uint i = 1; i < p2.size(); i++) p2[i] += p2[i-1];  // cumulative sum

    // We need to remove from `ints`
    while (d < 0) {
        long double u = static_cast<long double>(eng()) / samplers::pcg_max64;
        uint ind = std::lower_bound(p2.begin(), p2.end(), u) - p2.begin();
        ints[ind]--;
        d++;
    }
    // We need to add to `ints`
    while (d > 0) {
        long double u = static_cast<long double>(eng()) / samplers::pcg_max64;
        uint ind = std::lower_bound(p2.begin(), p2.end(), u) - p2.begin();
        ints[ind]++;
        d--;
    }

    return;
}


template<typename RNG>
class TableTable64 {
public:
    // Stores vectors of each category's Pr(sampled):
    std::vector<std::vector<uint>> T;
    // Stores values at which to transition between vectors of `T`:
    std::vector<uint64> t;

    TableTable64(const std::vector<long double>& probs, RNG& eng) : T(4), t(3, 0) {

        uint n_tables = T.size();

        uint n = probs.size();
        std::vector<uint64> ints(n);
        // Filling the `ints` vector based on `probs`
        fill_ints64<RNG>(probs, ints, eng);

        std::vector<uint> sizes(n_tables, 0);
        // Adding up sizes of `T` vectors:
        for (uint i = 0; i < n; i++) {
            for (uint k = 1; k <= n_tables; k++) {
                sizes[k-1] += dg64(ints[i], k);
            }
        }

        // Adding up thresholds in the `t` vector
        for (uint64 k = 0; k < (n_tables - 1); k++) {
            t[k] = static_cast<uint64>(sizes[k]);
            uint64 kk = 64 - 16 * (1 + k);
            t[k]<<=kk;
            if (k > 0) t[k] += t[k-1];
        }
        // Re-sizing `T` vectors:
        for (uint i = 0; i < n_tables; i++) T[i].resize(sizes[i]);

        // Filling `T` vectors
        for (uint k = 1; k <= n_tables; k++) {
            uint ind = 0; // index inside `T[k-1]`
            for (uint i = 0; i < n; i++) {
                uint z = dg64(ints[i], k);
                for (uint j = 0; j < z; j++) T[k-1][ind + j] = i;
                ind += z;
            }
        }
    }

    uint sample(RNG& eng) {
        uint64 j = eng();
        if (j<t[0]) return T[0][j>>(64-16*1)];
        if (j<t[1]) return T[1][(j-t[0])>>(64-16*2)];
        if (j<t[2]) return T[2][(j-t[1])>>(64-16*3)];
        return T[3][j-t[2]];
    }

private:
    static uint dg64(const uint64& m, const uint& k) {
        uint64 x = ((m>>(64ULL-16ULL*static_cast<uint64>(k)))&65535ULL);
        uint y = static_cast<uint>(x);
        return y;
    }
};


//[[Rcpp::export]]
SEXP make_table64(const std::vector<long double>& probs) {

    pcg64 eng = seeded_pcg64<pcg64>();

    XPtr<TableTable64<pcg64>> tt(new TableTable64<pcg64>(probs, eng), true);

    return tt;

}

//[[Rcpp::export]]
SEXP make_table64_fast(const std::vector<long double>& probs) {

    pcg64_fast eng = seeded_pcg64<pcg64_fast>();

    XPtr<TableTable64<pcg64_fast>> tt(new TableTable64<pcg64_fast>(probs, eng), true);

    return tt;
}


template <typename RNG>
uint sample_table_rare64_(SEXP xptr_sexp, const uint64& N, const uint& rare) {

    XPtr<TableTable64<RNG>> xptr(xptr_sexp);

    uint rares = 0;

    RNG eng = seeded_pcg64<RNG>();

    for (uint64 i = 0; i < N; i++) {
        uint k = xptr->sample(eng);
        if (k == rare) rares++;
    }

    return rares;
}


//[[Rcpp::export]]
uint sample_table_rare64(SEXP tt_, const uint64& N, const uint& rare) {

    uint rares = sample_table_rare64_<pcg64>(tt_, N, rare);

    return rares;
}

//[[Rcpp::export]]
uint sample_table_rare64_fast(SEXP tt_, const uint64& N, const uint& rare) {

    uint rares = sample_table_rare64_<pcg64_fast>(tt_, N, rare);

    return rares;
}





/*
 Fast Generation of Discrete Random Variables by Marsaglia, Tsang, & Wang (2004)
 JSS Journal of Statistical Software
 */


//[[Rcpp::export]]
SEXP make_table(const std::vector<double>& probs) {

    XPtr<TableTable> tt(new TableTable(probs), true);

    return tt;

}

//[[Rcpp::export]]
std::vector<uint> sample_table(SEXP tt_, const uint64& N) {
    XPtr<TableTable> tt(tt_);

    std::vector<uint> out(N);
    uint max_j = 0;

    pcg32 eng(static_cast<uint>(R::runif(0,1) * 2147483647));

    for (uint64 i = 0; i < N; i++) {
        out[i] = tt->sample(eng);
    }

    Rcout << max_j << std::endl;

    return out;
}


//[[Rcpp::export]]
uint sample_table_rare(SEXP tt_, const uint64& N, const uint& rare) {

    uint rares = sample_rare_<TableTable>(tt_, N, rare);

    return rares;
}

//[[Rcpp::export]]
List see_table(SEXP tt_) {
    XPtr<TableTable> tt(tt_);
    List out = List::create(_["AA"] = tt->T[0],
                            _["BB"] = tt->T[1],
                            _["CC"] = tt->T[2],
                            _["DD"] = tt->T[3],
                            _["t"] = tt->t);
    return out;
}





// alias sampling of unsigned integers
class AliasUInts {
public:
    AliasUInts() : F(), L(), n(0) {};
    AliasUInts(const std::vector<double>& p, const double& tol = SMALL_TOLERANCE);
    // To get the length of F (and L bc they should always be the same)
    uint size() const noexcept {
        return n;
    }
    // Actual alias sampling
    inline uint sample(pcg32& eng) const {
        // uniform in range [0,1)
        double u = static_cast<double>(eng()) / samplers::pcg_max;
        // Not doing +1 [as is done in Yang (2006)] to keep it in 0-based indexing
        uint k = n * u;
        double r = n * u - k;
        if (r >= F[k]) k = L[k];
        return k;
    };
private:
    std::vector<double> F;
    std::vector<uint> L;
    uint n;
};

AliasUInts::AliasUInts(const std::vector<double>& probs, const double& tol) {

    n = probs.size();

    // F_ and L_ are temporary vectors in arma formats to make the math easier
    arma::vec p(probs);
    double sum_p = arma::accu(p);
    p /= sum_p;
    arma::vec F_ = n * p;
    arma::uvec L_ = arma::regspace<arma::uvec>(0, n - 1);
    arma::ivec I(n);
    for (uint i = 0; i < n; i++) {
        if (F_(i) == 1) {
            L_(i) = i;
            I(i) = 0;
        } else if (F_(i) < 1) {
            I(i) = -1;
        } else {
            I(i) = 1;
        }
    }

    while (arma::any(I != 0)) {

        arma::uvec jv = arma::find(I == -1);  // underfull (i.e., F_ < 1)
        arma::uvec kv = arma::find(I == 1);  // overfull (i.e. F_ > 1)
        uint j = jv(0);
        if (kv.n_elem == 0) {
            stop("Numerical issue. Difference between one of the entries ",
                 "and 1 is " + std::to_string(F_(j) - 1));
        }
        uint k = kv(0);
        L_(j) = k;
        F_(k) = F_(k) - (1 - F_(j));
        I(j) = 0;
        if (std::abs(1 - F_(k)) < tol) {
            F_(k) = 1;
            I(k) = 0;
        } else if (F_(k) < 1) {
            I(k) = -1;
        }
    }

    F = arma::conv_to<std::vector<double>>::from(F_);
    L = arma::conv_to<std::vector<uint>>::from(L_);

    return;

}




class IndexTable {
public:

    IndexTable(const std::vector<double>& p_)
        : P(p_.size()), Q(p_.size()), m(p_.size()) {

        double sum_p = std::accumulate(p_.begin(), p_.end(), 0.0);

        std::vector<double> p(p_);

        for (uint i = 0; i < m; i++) {
            p[i] /= sum_p;
            P[i] = p[i];
            if (i > 0) P[i] += P[i-1];
        }
        uint i = 0;
        for (uint j = 0; j < m; j++) {
            while (P[i] < static_cast<double>(j) / static_cast<double>(m)) i++;
            Q[j] = i;
        }
    };

    uint sample(pcg32& eng) {
        double U = static_cast<double>(eng()) / samplers::pcg_max;
        uint j = m * U;
        uint i = Q[j];
        while (U >= P[i]) i++;
        return i;
    }

    void print() const {
        arma::rowvec x(P);
        x.print("P:");
        arma::urowvec y(Q);
        y.print("Q:");
        return;
    }

private:
    std::vector<double> P;
    std::vector<uint> Q;
    uint m;
};


//[[Rcpp::export]]
SEXP make_index(const std::vector<double>& p) {
    XPtr<IndexTable> it(new IndexTable(p), true);
    return it;
}

//[[Rcpp::export]]
std::vector<uint> sample_index(SEXP it_, const uint64& N) {
    XPtr<IndexTable> it(it_);

    std::vector<uint> out(N);

    pcg32 eng(static_cast<uint>(R::runif(0,1) * 2147483647));

    for (uint64 i = 0; i < N; i++) {
        out[i] = it->sample(eng);
    }
    return out;
}

//[[Rcpp::export]]
uint sample_index_rare(SEXP it_, const uint64& N, const uint& rare) {

    uint rares = sample_rare_<IndexTable>(it_, N, rare);

    return rares;
}

//[[Rcpp::export]]
void see_index(SEXP it_) {
    XPtr<IndexTable> it(it_);
    it->print();
    return;
}



//[[Rcpp::export]]
SEXP make_alias(const std::vector<double>& p, const double& tol) {
    XPtr<AliasUInts> aup(new AliasUInts(p, tol), true);
    return aup;
}

//[[Rcpp::export]]
std::vector<uint> sample_alias(SEXP aup_, const uint64& N) {
    XPtr<AliasUInts> aup(aup_);

    std::vector<uint> out(N);

    pcg32 eng(static_cast<uint>(R::runif(0,1) * 2147483647));

    for (uint64 i = 0; i < N; i++) {
        out[i] = aup->sample(eng);
    }
    return out;
}

//[[Rcpp::export]]
uint sample_alias_rare(SEXP aup_, const uint64& N, const uint& rare) {

    uint rares = sample_rare_<AliasUInts>(aup_, N, rare);

    return rares;
}


