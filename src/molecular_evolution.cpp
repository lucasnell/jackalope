#include <RcppArmadillo.h>
#include <sitmo.h>    // sitmo
#include <vector>  // vector class
#include <string>  // string class
#include <algorithm>  // lower_bound, random_shuffle
#include <cmath>  // std::exp, std::log
#include <numeric>  // accumulate
#include <unordered_map>  // unordered_map
#include <queue>  // priority_queue
#include <deque>  // deque
#ifdef _OPENMP
#include <omp.h>  // omp
#endif


#include "gemino_types.h"
#include "molecular_evolution.h"
#include "alias.h"    // Alias sampling

using namespace Rcpp;




// Get a sequence's overall event rate
//[[Rcpp::export]]
double seq_rate(std::string seq_, const XPtr<QMaps> qm) {
    double rate = 0;
    for (char c : seq_) {
        rate += qm->q[c];
    }
    return rate;
}





// Generate a string to test sampling characters from a string
//[[Rcpp::export]]
std::string gen_test_str(uint len) {

    std::string out = "";
    uint n_each = len / 4;

    for (char c : mevo::bases) {
        for (uint i = 0; i < n_each; i++) {
            out += c;
        }
    }

    // std::random_shuffle(out.begin(), out.end());

    return out;
}



// equivalent to a^b (this way is faster than using std::pow(a, b))
inline double fast_pow(double a, double b) {
    return std::exp(b * std::log(a));
}

// uniform in range (0,1)
inline double runif_01(sitmo::prng_engine& eng) {
    return (static_cast<double>(eng()) + 1) / (mevo::sitmo_max + 1);
}
// uniform in range (a,b)
inline double runif_ab(sitmo::prng_engine& eng, double a, double b) {
    return a + (static_cast<double>(eng()) / (mevo::sitmo_max + 1)) * (b - a);
}


/*
 At a time when an event occurs, this samples nucleotides based on their rates
 and returns a random location where the event will occur

 Reservoir sampling used from...
 Efraimidis, P. S., and P. G. Spirakis. 2006. Weighted random sampling with a
 reservoir. Information Processing Letters 97:181–185.
 */
//[[Rcpp::export]]
uint event_location(const std::string& S,
                    XPtr<QMaps> qm, uint seed) {

    if (S.size() == 0) stop("Empty string sent to event_location.");
    if (S.size() == 1) return 0;

    sitmo::prng_engine eng(seed);

    double r, key, X, w, t;
    uint N = S.size();

    // Create a NucleoKeys object to store position and key
    r = runif_01(eng);
    key = std::log(r) / qm->w[S[0]];
    key = std::exp(key);
    NucleoKeyPos pq(key, 0);

    uint c = 0;
    while (c < (N-1)) {
        r = runif_01(eng);
        X = std::log(r) / std::log(pq.key);
        uint i = c + 1;
        double wt_sum0 = qm->w[S[c]];
        double wt_sum1 = qm->w[S[c]] + qm->w[S[i]];
        while (X > wt_sum1 && i < (N-1)) {
            i++;
            wt_sum0 += qm->w[S[(i-1)]];
            wt_sum1 += qm->w[S[i]];
        }
        if (X > wt_sum1) break;
        if (wt_sum0 >= X) continue;

        pq.pos = i;

        w = qm->w[S[i]];
        t = fast_pow(pq.key, w);
        r = runif_ab(eng, t, 1.0);
        key = fast_pow(r, 1 / w);
        pq.key = key;

        c = i;
    }

    return pq.pos;
}




/*
 Get weight when you have a reference sequence and variant sites/nucleos
 This version also takes care of "iterating": `++` for indices and `.pop_front()` for
 deque objects.
 NEVER TESTED --> NOT GUARANTEED TO WORK
 */
double variant_event_wt(const std::string& S,
                        std::deque<uint>& sites,
                        std::deque<std::string>& nucleos,
                        uint& ind, uint& ni,
                        const XPtr<QMaps> qm) {

    double wt;

    if (sites.empty()) {
        wt = qm->w[S[ind]];
        ind++;
    } else if (sites.front() == ind) {
        while (nucleos.front().size() == 0 && sites.front() == ind) {
            sites.pop_front();
            nucleos.pop_front();
            ind++;
        }
        if (sites.front() != ind) {
            wt = qm->w[S[ind]];
            ind++;
        } else {
            wt = qm->w[nucleos.front()[ni]];
            if (ni == (nucleos.front().size() - 1)) {
                sites.pop_front();
                nucleos.pop_front();
                ni = 0;
                ind++;
            } else {
                ni++;
            }
        }
    } else {
        wt = qm->w[S[ind]];
        ind++;
    }
    return wt;
}

/*
 This is another version of above for when you don't want indices or deque objects
 to be "iterated": `++` for indices and `pop_front` for deque objects.
 NEVER TESTED --> NOT GUARANTEED TO WORK
 */
double variant_event_wt_noit(const std::string& S,
                             const std::deque<uint>& sites,
                             const std::deque<std::string>& nucleos,
                             const uint& ind, const uint& ni,
                             const XPtr<QMaps> qm) {

    double wt;
    // Internal copies of indices that can be iterated:
    uint indi = ind, nii = ni;
    // Index for nucleos and sites objects to use instead of `.front()`
    uint ii = 0;

    if (sites.empty()) {
        wt = qm->w[S[ind]];
    } else if (sites[ii] == ind) {
        while (nucleos[ii].size() == 0 && sites[ii] == ind) {
            ii++;
            indi++;
            if (ii >= sites.size()) {
                wt = qm->w[S[ind]];
            }
        }
        if (sites[ii] != ind) {
            wt = qm->w[S[ind]];
            indi++;
        } else {
            wt = qm->w[nucleos[ii][ni]];
            if (ni == (nucleos[ii].size() - 1)) {
                ii++;
                nii = 0;
                indi++;
            } else {
                nii++;
            }
        }
    } else {
        wt = qm->w[S[ind]];
        indi++;
    }
    return wt;
}


/*
 This is the same as above, but uses variant info
 */
//[[Rcpp::export]]
uint event_location2(const std::string& S,
                     std::deque<uint> sites,
                     std::deque<std::string> nucleos,
                     XPtr<QMaps> qm, uint seed) {

    if (S.size() == 0 && nucleos.size() == 0) {
        stop("Empty string sent to event_location.");
    }
    if (S.size() == 1) return 0;

    if (nucleos.size() != sites.size()) stop("sites and nucleos aren't same length");

    sitmo::prng_engine eng(seed);

    double r, key, X, w, t;
    uint N = S.size();

    // `c` is for the position in `S`
    // `ni` is for the position within a single, multi-character `nucleos` string
    uint c = 0, ni = 0;
    // Create a NucleoKeys object to store position and key
    w = variant_event_wt(S, sites, nucleos, c, ni, qm);
    r = runif_01(eng);
    key = std::log(r) / w;
    key = std::exp(key);
    NucleoKeyPos pq(key, c);

    while (c < (N-1)) {
        r = runif_01(eng);
        X = std::log(r) / std::log(pq.key);
        uint i = c + 1;
        double wt_sum0 = variant_event_wt(S, sites, nucleos, c, ni, qm);
        double wt_sum1 = wt_sum0 + variant_event_wt(S, sites, nucleos, i, ni, qm);
        while (X > wt_sum1 && i < (N-1)) {
            i++;
            // LEFT OFF BELOW --> HOW TO GET i-1 WORKING WITH FXN ABOVE??
            wt_sum0 += variant_event_wt(S, sites, nucleos, i-1, ni, qm);
            wt_sum1 += qm->w[S[i]];
        }
        if (X > wt_sum1) break;
        if (wt_sum0 >= X) continue;

        pq.pos = i;

        w = qm->w[S[i]];
        t = fast_pow(pq.key, w);
        r = runif_ab(eng, t, 1.0);
        key = fast_pow(r, 1 / w);
        pq.key = key;

        c = i;
    }

    return pq.pos;
}


//[[Rcpp::export]]
arma::mat qm_Q(XPtr<QMaps> qm) {

    arma::mat Q(4, qm->Q['A'].size());
    uint i = 0;
    for (char c : mevo::bases) {
        arma::rowvec V(qm->Q[c]);
        Q.row(i) = V;
        i++;
    }
    return Q;
}

//[[Rcpp::export]]
arma::mat qm_M(XPtr<QMaps> qm) {

    arma::mat M(4, qm->M['A'].size());
    uint i = 0;
    for (char c : mevo::bases) {
        arma::rowvec V(qm->M[c]);
        M.row(i) = V;
        i++;
    }
    return M;
}

//[[Rcpp::export]]
std::vector<double> qm_q(XPtr<QMaps> qm) {

    std::vector<double> q(4);
    uint i = 0;
    for (char c : mevo::bases) {
        q[i] = qm->q[c];
        i++;
    }
    return q;
}



/*
 Q and M matrices, for rates and Pr(event), respectively, for a given nucleotide
 Also q, which is for nucleotides' rates
 They are stored as unordered_map objects for speed
 */
//[[Rcpp::export]]
XPtr<QMaps> QM(
        const double& beta, const double& alpha_1, const double& alpha_2,
        const double& xi, const double& psi,
        const double& pi_a, const double& pi_c,
        const double& pi_g, const double& pi_t) {

    // Number of event types: 4 nucleotide types, plus insertions and deletions from
    // 1 to 10 bp long
    uint n_events = 4 + 2 * 10;


    XPtr<QMaps> qm(new QMaps(n_events), true);


    double pi_y = pi_t + pi_c;
    double pi_r = pi_a + pi_g;

    // -----------
    // Rates
    // -----------
    qm->Q['A'][0] = -(alpha_2 * pi_g + beta * pi_y + xi);
    qm->Q['A'][1] = beta * pi_c;
    qm->Q['A'][2] = alpha_2 * pi_g;
    qm->Q['A'][3] = beta * pi_t;

    qm->Q['C'][0] = beta * pi_a;
    qm->Q['C'][1] = -(alpha_1 * pi_t + beta * pi_r + xi);
    qm->Q['C'][2] = beta * pi_g;
    qm->Q['C'][3] = alpha_1 * pi_t;

    qm->Q['G'][0] = alpha_2 * pi_a;
    qm->Q['G'][1] = beta * pi_c;
    qm->Q['G'][2] = -(alpha_2 * pi_a + beta * pi_y + xi);
    qm->Q['G'][3] = beta * pi_t;

    qm->Q['T'][0] = beta * pi_a;
    qm->Q['T'][1] = alpha_1 * pi_c;
    qm->Q['T'][2] = beta * pi_g;
    qm->Q['T'][3] = -(alpha_1 * pi_c + beta * pi_r + xi);

    // Adding columns for indels
    double xi_i = xi / (1 + 1/psi);  // insertion rate
    double xi_d = xi / (1 + psi);    // deletion rate
    // rates for indel lengths 1:10, that sum to 1
    arma::rowvec omega = arma::regspace<arma::rowvec>(-1, -1, -10);
    omega = arma::exp(omega);
    omega /= arma::accu(omega);
    for (char c : mevo::bases) {
        for (uint i = 0; i < omega.n_elem; i++) {
            qm->Q[c][i + 4] = omega(i) * xi_i;
            qm->Q[c][i + 4 + 10] = omega(i) * xi_d;
        }
    }

    // -----------
    // Event matrix (M_i = Pr(event_i))
    // and
    // Rates by nucleotide (q_i = rate for nucleotide_i)
    // -----------
    double q_sum = 0;
    for (char c : mevo::bases) {
        qm->M[c] = std::vector<double>(n_events);
        uint ind = mevo::base_inds[c];
        double qi = -(qm->Q[c][ind]);
        for (uint i = 0; i < n_events; i++) qm->M[c][i] = qm->Q[c][i] / qi;
        qm->M[c][ind] = 0;
        // rate for nucleotide `c`:
        qm->q[c] = qi;
        q_sum += qi;
    }
    // so q1 rates sum to one (for weighted sampling):
    for (char c : mevo::bases) qm->w[c] = qm->q[c] / q_sum;

    return qm;

}


/*
 Create the L and F vectors for alias sampling.

 For more info, see p 299 in...

 Yang, Z. 2006. Computational molecular evolution. (P. H. Harvey and R. M. May, Eds.).
     Oxford University Press, New York, NY, USA.
 */

//[[Rcpp::export]]
std::vector<std::vector<double>> alias_FL(arma::vec p,
                                          double tol = 0.00000001490116119385) {

    uint n = p.n_elem;

    arma::vec F = n * p;
    arma::ivec L = arma::regspace<arma::ivec>(0, n - 1);;
    arma::ivec I(n);
    for (uint i = 0; i < n; i++) {
        if (F(i) == 1) {
            L(i) = i;
            I(i) = 0;
        } else if (F(i) < 1) {
            I(i) = -1;
        } else {
            I(i) = 1;
        }
    }

    while (arma::any(I != 0)) {

        arma::uvec jv = arma::find(I == -1);  // underfull (i.e., F < 1)
        arma::uvec kv = arma::find(I == 1);  // overfull (i.e. F > 1)
        uint j = jv(0);
        if (kv.n_elem == 0) {
            stop("Numerical issue. Difference between one of the entries ",
                 "and 1 is "); // , F(j) - 1);
        }
        uint k = kv(0);
        L(j) = k;
        F(k) = F(k) - (1 - F(j));
        I(j) = 0;
        if (std::abs(1 - F(k)) < tol) {
            F(k) = 1;
            I(k) = 0;
        } else if (F(k) < 1) {
            I(k) = -1;
        }
    }

    std::vector<std::vector<double>> out(2);
    out[0] = arma::conv_to<std::vector<double>>::from(F);
    out[1] = arma::conv_to<std::vector<double>>::from(L);

    return out;

}




//[[Rcpp::export]]
XPtr<Alias> alias_tables(XPtr<QMaps> qm, std::vector<double> pis,
                         arma::umat event_matrix,
                         double tol = 0.00000001490116119385) {

    EventAliasTables am;
    am.event_types = arma::conv_to<std::vector<uint>>::from(event_matrix.row(0));
    am.event_lengths = arma::conv_to<std::vector<uint>>::from(event_matrix.row(1));

    // For sampling events that depend on the current nucleotide:
    for (char c : mevo::bases) {

        arma::vec p_c(qm->M[c]);

        std::vector<std::vector<double>> FL = alias_FL(p_c, tol);

        am.F[c] = FL[0];
        std::vector<uint> L(FL[1].begin(), FL[1].end());
        am.L[c] = L;
    }

    // For sampling new sequences (for insertions) that don't depend on a
    // current nucleotide:
    NucleoAliasTables amp;
    if (pis.size() != 4) stop("pi vector should be length 4");
    std::vector<std::vector<double>> FL = alias_FL(pis, tol);
    amp.F = FL[0];
    std::vector<uint> L(FL[1].begin(), FL[1].end());
    amp.L = L;


    XPtr<Alias> al(new Alias(am, amp), true);

    return al;

}


//[[Rcpp::export]]
arma::umat am_L(XPtr<Alias> al) {

    EventAliasTables& am(al->event);

    arma::umat M(4, am.L['A'].size());
    uint i = 0;
    for (char c : mevo::bases) {
        arma::urowvec V(am.L[c]);
        M.row(i) = V;
        i++;
    }
    return M;
}
//[[Rcpp::export]]
arma::mat am_F(XPtr<Alias> al) {

    EventAliasTables& am(al->event);

    arma::mat M(4, am.F['A'].size());
    uint i = 0;
    for (char c : mevo::bases) {
        arma::rowvec V(am.F[c]);
        M.row(i) = V;
        i++;
    }
    return M;
}


// uniform in range [0,1)
inline double runif_001(sitmo::prng_engine& eng) {
    return ((double) eng()) / (mevo::sitmo_max + 1);
}

/*
 Alias sampling for which event to perform
 */
//[[Rcpp::export]]
uint event_samp(char c, XPtr<Alias> al, uint seed) {

    sitmo::prng_engine eng(seed);

    double n = static_cast<double>(al->event.F[c].size());

    double u = runif_001(eng);
    // Not doing +1 [as is done in Yang (2006)] to keep it in 0-based indexing
    double dbl_k = std::floor(n * u);
    double r = n * u - dbl_k;
    uint k = static_cast<uint>(dbl_k);

    if (r < al->event.F[c][k]) return k;

    return al->event.L[c][k];

}


/*
 Alias sampling for which sequence to insert
*/
//[[Rcpp::export]]
std::string nucleo_samp(uint N, XPtr<Alias> al, uint seed) {

    sitmo::prng_engine eng(seed);

    std::string out(N, 'x');

    double n = 4;

    for (uint i = 0; i < N; i++) {
        double u = runif_001(eng);
        // Not doing +1 [as is done in Yang (2006)] to keep it in 0-based indexing
        double dbl_k = std::floor(n * u);
        double r = n * u - dbl_k;
        uint k = static_cast<uint>(dbl_k);
        if (r >= al->nucleo.F[k]) k = al->nucleo.L[k];
        out[i] = mevo::bases[k];
    }

    return out;
}







//' Consolidate reference and variant info into a single string.
//'
//' This should replace the current version once the new variant structure is
//' implemented.
//'
//' This function assumes nucleos and sites are both sorted by sites
//' (increasing order).
//'
//' @param str
//' @param sites
//' @param nucleos
//'
//' @return
//'
//' @noRd
//'
//'
//' @examples
//'
//'
//[[Rcpp::export]]
std::string consolidate(const std::string& str,
                        std::deque<uint> sites,
                        std::deque<std::string> nucleos) {
    std::string out = "";
    uint i = 0;
    while (!sites.empty()) {
        for (uint j = i; j < sites.front(); j++) {
            out += str[j];
        }
        i = sites.front();
        while (i == sites.front()) {
            out += nucleos.front();
            nucleos.pop_front();
            sites.pop_front();
            i++;
            if (sites.empty()) break;
        }
    }
    for (uint j = i; j < str.size(); j++) out += str[j];

    return out;
}






//' Combine two sets of variant info into a single set.
//'
//' This should replace the current version once the new variant structure is
//' implemented.
//'
//' This function assumes nucleos and sites are both sorted by sites
//' (increasing order).
//'
//' @param sites_orig
//' @param nucleos_orig
//' @param sites_new
//' @param nucleos_new
//' @param qm
//' @param seed
//'
//' @return
//'
//' @noRd
//'
//'
//' @examples
//'
//'
//[[Rcpp::export]]
List combine_variants(std::deque<uint> sites_orig,
                      std::deque<std::string> nucleos_orig,
                      std::deque<uint> sites_new,
                      std::deque<std::string> nucleos_new,
                      const XPtr<QMaps> qm, const uint& seed) {

    if (sites_orig.size() != nucleos_orig.size()) {
        stop("Original sites and nucleos vectors are not the same length.");
    } else if (sites_new.size() != nucleos_new.size()) {
        stop("New sites and nucleos vectors are not the same length.");
    }

    std::deque<uint> sites;
    std::deque<std::string> nucleos;

    while (!sites_orig.empty() || !sites_new.empty()) {
        Rcpp::checkUserInterrupt();
        if (sites_orig.empty()) {
            sites.push_back(sites_new.front());
            nucleos.push_back(nucleos_new.front());
            sites_new.pop_front();
            nucleos_new.pop_front();
        } else if (sites_new.empty()) {
            sites.push_back(sites_orig.front());
            nucleos.push_back(nucleos_orig.front());
            sites_orig.pop_front();
            nucleos_orig.pop_front();
        } else if (sites_new.front() > sites_orig.front()) {
            sites.push_back(sites_orig.front());
            nucleos.push_back(nucleos_orig.front());
            sites_orig.pop_front();
            nucleos_orig.pop_front();
        } else if (sites_new.front() < sites_orig.front()) {
            sites.push_back(sites_new.front());
            nucleos.push_back(nucleos_new.front());
            sites_new.pop_front();
            nucleos_new.pop_front();
        // If `sites_new[i]` is in `sites_orig`, then we have to add it appropriately
        } else {
            std::string combined_str;
            if (nucleos_orig.front().size() == 0 ||
                (nucleos_new.front().size() <= 1 &&
                nucleos_orig.front().size() == 1)) {
                combined_str = nucleos_new.front();
            } else {
                // Sampling of locations within the original insertion.
                // (If `nucleos_orig.front().size() == 1`, `j` will always be 0.)
                // This will be the location for the new mutation
                uint j = event_location(nucleos_orig.front(), qm, seed);

                std::string& nt0(nucleos_orig.front());
                std::string& nt1(nucleos_new.front());
                combined_str = "";
                // Bc insertions includes the original reference sequence,
                // we want to replace that with the original *insertion* sequence now.
                if (nt1.size() > 1) nt1[0] = nt0[j];
                // Combining the original and new sequences now:
                for (uint k = 0; k < j; k++) combined_str += nt0[k];
                combined_str += nt1;
                for (uint k = j+1; k < nt0.size(); k++) {
                    combined_str += nt0[k];
                }
            }

            nucleos.push_back(combined_str);
            sites.push_back(sites_new.front());

            sites_new.pop_front();
            nucleos_new.pop_front();
            sites_orig.pop_front();
            nucleos_orig.pop_front();
        }
    }

    return List::create(_["sites"] = sites, _["nucleos"] = nucleos);
}
