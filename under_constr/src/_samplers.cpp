#include <RcppArmadillo.h>
#include <sitmo.h>    // sitmo prng
#include <vector>
#include <string>



//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(sitmo, RcppArmadillo, gemino)]]

using namespace Rcpp;


#define SMALL_TOLERANCE 0.000000000001
#define MAX_UINT 4294967295

typedef uint_fast32_t uint;

namespace test {
    double sitmo_max = static_cast<double>(sitmo::prng_engine::max());
}




uint dg(uint m, uint k) {
    // uint x = ((m>>(30-6*k))&63);
    uint x = ((m>>(32-4*k))&15);
    return x;
}


//[[Rcpp::export]]
uint dg_R(uint m, uint k) {
    uint x = dg(m, k);
    return x;
}


class Table {
public:
    std::vector<uint> AA;
    std::vector<uint> BB;
    std::vector<uint> CC;
    std::vector<uint> DD;
    std::vector<uint> EE;
    std::vector<uint> FF;
    std::vector<uint> GG;
    std::vector<uint> HH;

    std::vector<uint>& operator[](const uint& idx) {
        switch(idx) {
        case 0 : return AA;
        case 1 : return BB;
        case 2 : return CC;
        case 3 : return DD;
        case 4 : return EE;
        case 5 : return FF;
        case 6 : return GG;
        default : return HH;
        }
    }
    const std::vector<uint>& operator[](const uint& idx) const {
        switch(idx) {
        case 0 : return AA;
        case 1 : return BB;
        case 2 : return CC;
        case 3 : return DD;
        case 4 : return EE;
        case 5 : return FF;
        case 6 : return GG;
        default : return HH;
        }
    }

    uint size() const {
        return size_;
    }

private:
    uint size_ = 8;
};


struct TableTable {
    Table T;
    std::vector<uint> t;

    TableTable(const std::vector<double>& probs) : T(), t(7, 0) {

        uint n_tables = T.size();

        uint n = probs.size();
        std::vector<uint> ints(n);
        std::vector<uint> sizes(n_tables, 0);
        for (uint i = 0; i < n; i++) {
            double tmp = std::round(probs[i] * 4294967296);
            if (tmp > MAX_UINT) tmp = MAX_UINT;
            ints[i] = static_cast<uint>(tmp);
            for (uint k = 1; k <= n_tables; k++) {
                sizes[k-1] += dg(ints[i], k);
            }
        }
        for (uint k = 0; k < (n_tables - 1); k++) {
            // t[k] = sizes[k]<<(30-6*(1+k));
            t[k] = sizes[k]<<(32-4*(1+k));
            if (k > 0) t[k] += t[k-1];
        }
        for (uint i = 0; i < n_tables; i++) T[i].resize(sizes[i]);
        std::vector<uint> inds(n_tables, 0);
        for (uint i = 0; i < n; i++){
            uint& m = ints[i];
            for (uint k = 1; k <= n_tables; k++) {
                uint z = dg(m,k);
                for (uint j = 0; j < z; j++) T[k-1][inds[k-1] + j] = i;
                inds[k-1] += z;
            }
        }
    }

    uint sample(sitmo::prng_engine& eng, uint& max_j) {
        uint j = eng();
        // j >>= 2; // going from 32-bit to 30-bit
        // j /= 4; // going from 32-bit to 30-bit
        if (j > max_j) max_j = j;
        if (j<t[0]) return T[0][j>>28];
        for (uint i = 1; i < (T.size() - 1); i++) {
            if (j<t[i]) return T[i][(j-t[i-1])>>(32-4*(i+1))];
        }
        return T[T.size()-1][j-t[T.size()-2]];
        // if (j<t[0]) return T.AA[j>>24];
        // if (j<t[1]) return T.BB[(j-t[0])>>18];
        // if (j<t[2]) return T.CC[(j-t[1])>>12];
        // if (j<t[3]) return T.DD[(j-t[2])>>6];
        // return T.EE[j-t[3]];
    }

    void print() const {
        std::vector<std::string> names = {"AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH"};
        for (uint i = 0; i < T.size(); i++) {
            arma::urowvec x(T[i]);
            x.print(names[i] + ":");
        }
        arma::urowvec x(t);
        x.print("t:");
    }

};




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
std::vector<uint> sample_table(SEXP tt_, const uint& N) {
    XPtr<TableTable> tt(tt_);

    std::vector<uint> out(N);
    uint max_j = 0;

    sitmo::prng_engine eng(static_cast<uint>(R::runif(0,1) * 2147483647));

    for (uint i = 0; i < N; i++) {
        out[i] = tt->sample(eng, max_j);
    }

    Rcout << max_j << std::endl;

    return out;
}


//[[Rcpp::export]]
uint sample_table_rare(SEXP tt_, const uint& N, const uint& rare) {
    XPtr<TableTable> tt(tt_);

    uint rares = 0;
    uint max_j = 0;

    sitmo::prng_engine eng(static_cast<uint>(R::runif(0,1) * 2147483647));

    for (uint i = 0; i < N; i++) {
        uint k = tt->sample(eng, max_j);
        if (k == rare) rares++;
    }
    Rcout << max_j << std::endl;
    return rares;
}

//[[Rcpp::export]]
void see_table(SEXP tt_) {
    XPtr<TableTable> tt(tt_);
    tt->print();
    return;
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
    inline uint sample(sitmo::prng_engine& eng) const {
        // uniform in range [0,1)
        double u = static_cast<double>(eng()) / test::sitmo_max;
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

    uint sample(sitmo::prng_engine& eng) {
        double U = static_cast<double>(eng()) / test::sitmo_max;
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
std::vector<uint> sample_index(SEXP it_, const uint& N) {
    XPtr<IndexTable> it(it_);

    std::vector<uint> out(N);

    sitmo::prng_engine eng(static_cast<uint>(R::runif(0,1) * 2147483647));

    for (uint i = 0; i < N; i++) {
        out[i] = it->sample(eng);
    }
    return out;
}

//[[Rcpp::export]]
uint sample_index_rare(SEXP it_, const uint& N, const uint& rare) {
    XPtr<IndexTable> it(it_);

    uint rares = 0;

    sitmo::prng_engine eng(static_cast<uint>(R::runif(0,1) * 2147483647));

    for (uint i = 0; i < N; i++) {
        uint k = it->sample(eng);
        if (k == rare) rares++;
    }
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
std::vector<uint> sample_alias(SEXP aup_, const uint& N) {
    XPtr<AliasUInts> aup(aup_);

    std::vector<uint> out(N);

    sitmo::prng_engine eng(static_cast<uint>(R::runif(0,1) * 2147483647));

    for (uint i = 0; i < N; i++) {
        out[i] = aup->sample(eng);
    }
    return out;
}

//[[Rcpp::export]]
uint sample_alias_rare(SEXP aup_, const uint& N, const uint& rare) {
    uint rares = 0;
    XPtr<AliasUInts> aup(aup_);

    sitmo::prng_engine eng(static_cast<uint>(R::runif(0,1) * 2147483647));

    for (uint i = 0; i < N; i++) {
        uint k = aup->sample(eng);
        if (k == rare) rares++;
    }

    return rares;
}


