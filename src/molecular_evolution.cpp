#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <unordered_map>  // unordered_map
#include <deque>  // deque
#ifdef _OPENMP
#include <omp.h>  // omp
#endif


#include "gemino_types.h"
#include "molecular_evolution.h"
#include "sequence_classes.h"  // Var* and Ref* classes
#include "pcg.h"  // pcg seeding
#include "table_sampler.h"  // table method of sampling

using namespace Rcpp;




// // Get a sequence's overall event rate
// //[[Rcpp::export]]
// double seq_rate(std::string seq_, const XPtr<QMaps> qm) {
//     double rate = 0;
//     for (char c : seq_) {
//         rate += qm->q[c];
//     }
//     return rate;
// }
//
//
//
//
// /*
//  At a time when an event occurs, this samples nucleotides based on their rates
//  and returns a random location where the event will occur
//
//  Reservoir sampling used from...
//  Efraimidis, P. S., and P. G. Spirakis. 2006. Weighted random sampling with a
//  reservoir. Information Processing Letters 97:181–185.
//  */
// //[[Rcpp::export]]
// uint event_location(const std::string& S,
//                     XPtr<QMaps> qm, uint seed) {
//
//     if (S.size() == 0) stop("Empty string sent to event_location.");
//     if (S.size() == 1) return 0;
//
//     pcg32 eng(seed);
//
//     double r, key, X, w, t;
//     uint N = S.size();
//
//     // Create a NucleoKeys object to store position and key
//     r = runif_01(eng);
//     key = std::log(r) / qm->w[S[0]];
//     key = std::exp(key);
//     NucleoKeyPos pq(key, 0);
//
//     uint c = 0;
//     while (c < (N-1)) {
//         r = runif_01(eng);
//         X = std::log(r) / std::log(pq.key);
//         uint i = c + 1;
//         double wt_sum0 = qm->w[S[c]];
//         double wt_sum1 = qm->w[S[c]] + qm->w[S[i]];
//         while (X > wt_sum1 && i < (N-1)) {
//             i++;
//             wt_sum0 += qm->w[S[(i-1)]];
//             wt_sum1 += qm->w[S[i]];
//         }
//         if (X > wt_sum1) break;
//         if (wt_sum0 >= X) continue;
//
//         pq.pos = i;
//
//         w = qm->w[S[i]];
//         t = std::pow(pq.key, w);
//         r = runif_ab(eng, t, 1.0);
//         key = std::pow(r, 1 / w);
//         pq.key = key;
//
//         c = i;
//     }
//
//     return pq.pos;
// }
//
//
//
//
// /*
//  Get weight when you have a reference sequence and variant sites/nucleos
//  This version also takes care of "iterating": `++` for indices and `.pop_front()` for
//  deque objects.
//  NEVER TESTED --> NOT GUARANTEED TO WORK
//  */
// double variant_event_wt(const std::string& S,
//                         std::deque<uint>& sites,
//                         std::deque<std::string>& nucleos,
//                         uint& ind, uint& ni,
//                         const XPtr<QMaps> qm) {
//
//     double wt;
//
//     if (sites.empty()) {
//         wt = qm->w[S[ind]];
//         ind++;
//     } else if (sites.front() == ind) {
//         while (nucleos.front().size() == 0 && sites.front() == ind) {
//             sites.pop_front();
//             nucleos.pop_front();
//             ind++;
//         }
//         if (sites.front() != ind) {
//             wt = qm->w[S[ind]];
//             ind++;
//         } else {
//             wt = qm->w[nucleos.front()[ni]];
//             if (ni == (nucleos.front().size() - 1)) {
//                 sites.pop_front();
//                 nucleos.pop_front();
//                 ni = 0;
//                 ind++;
//             } else {
//                 ni++;
//             }
//         }
//     } else {
//         wt = qm->w[S[ind]];
//         ind++;
//     }
//     return wt;
// }
//
// /*
//  This is another version of above for when you don't want indices or deque objects
//  to be "iterated": `++` for indices and `pop_front` for deque objects.
//  NEVER TESTED --> NOT GUARANTEED TO WORK
//  */
// double variant_event_wt_noit(const std::string& S,
//                              const std::deque<uint>& sites,
//                              const std::deque<std::string>& nucleos,
//                              const uint& ind, const uint& ni,
//                              const XPtr<QMaps> qm) {
//
//     double wt;
//     // Internal copies of indices that can be iterated:
//     uint indi = ind, nii = ni;
//     // Index for nucleos and sites objects to use instead of `.front()`
//     uint ii = 0;
//
//     if (sites.empty()) {
//         wt = qm->w[S[ind]];
//     } else if (sites[ii] == ind) {
//         while (nucleos[ii].size() == 0 && sites[ii] == ind) {
//             ii++;
//             indi++;
//             if (ii >= sites.size()) {
//                 wt = qm->w[S[ind]];
//             }
//         }
//         if (sites[ii] != ind) {
//             wt = qm->w[S[ind]];
//             indi++;
//         } else {
//             wt = qm->w[nucleos[ii][ni]];
//             if (ni == (nucleos[ii].size() - 1)) {
//                 ii++;
//                 nii = 0;
//                 indi++;
//             } else {
//                 nii++;
//             }
//         }
//     } else {
//         wt = qm->w[S[ind]];
//         indi++;
//     }
//     return wt;
// }
//
//
// /*
//  This is the same as above, but uses variant info
//  */
// //[[Rcpp::export]]
// uint event_location2(const std::string& S,
//                      std::deque<uint> sites,
//                      std::deque<std::string> nucleos,
//                      XPtr<QMaps> qm, uint seed) {
//
//     if (S.size() == 0 && nucleos.size() == 0) {
//         stop("Empty string sent to event_location.");
//     }
//     if (S.size() == 1) return 0;
//
//     if (nucleos.size() != sites.size()) stop("sites and nucleos aren't same length");
//
//     pcg32 eng(seed);
//
//     double r, key, X, w, t;
//     uint N = S.size();
//
//     // `c` is for the position in `S`
//     // `ni` is for the position within a single, multi-character `nucleos` string
//     uint c = 0, ni = 0;
//     // Create a NucleoKeys object to store position and key
//     w = variant_event_wt(S, sites, nucleos, c, ni, qm);
//     r = runif_01(eng);
//     key = std::log(r) / w;
//     key = std::exp(key);
//     NucleoKeyPos pq(key, c);
//
//     while (c < (N-1)) {
//         r = runif_01(eng);
//         X = std::log(r) / std::log(pq.key);
//         uint i = c + 1;
//         double wt_sum0 = variant_event_wt(S, sites, nucleos, c, ni, qm);
//         double wt_sum1 = wt_sum0 + variant_event_wt(S, sites, nucleos, i, ni, qm);
//         while (X > wt_sum1 && i < (N-1)) {
//             i++;
//             // LEFT OFF BELOW --> HOW TO GET i-1 WORKING WITH FXN ABOVE??
//             wt_sum0 += variant_event_wt(S, sites, nucleos, i-1, ni, qm);
//             wt_sum1 += qm->w[S[i]];
//         }
//         if (X > wt_sum1) break;
//         if (wt_sum0 >= X) continue;
//
//         pq.pos = i;
//
//         w = qm->w[S[i]];
//         t = std::pow(pq.key, w);
//         r = runif_ab(eng, t, 1.0);
//         key = std::pow(r, 1 / w);
//         pq.key = key;
//
//         c = i;
//     }
//
//     return pq.pos;
// }





//' Q matrix for rates for a given nucleotide using the TN93 substitution model.
//'
//' @noRd
//'
std::unordered_map<char, std::vector<double>> TN93_rate_matrix(
        const double& pi_a, const double& pi_c,
        const double& pi_g, const double& pi_t,
        const double& alpha_1, const double& alpha_2, const double& beta,
        const double& xi) {

    double pi_y = pi_t + pi_c;
    double pi_r = pi_a + pi_g;

    /*
     Rate vectors ("Q" matrix in Yang (2006)).
     (These vectors are only used in the construction of the final object, so they don't
     need to be retained in their own class.)
    */
    std::unordered_map<char, std::vector<double>> Q;
    Q = {{'A', std::vector<double>(4, 0.0)},
         {'C', std::vector<double>(4, 0.0)},
         {'G', std::vector<double>(4, 0.0)},
         {'T', std::vector<double>(4, 0.0)}};
    Q['A'][0] = -(alpha_2 * pi_g + beta * pi_y + xi);
    Q['A'][1] = beta * pi_c;
    Q['A'][2] = alpha_2 * pi_g;
    Q['A'][3] = beta * pi_t;

    Q['C'][0] = beta * pi_a;
    Q['C'][1] = -(alpha_1 * pi_t + beta * pi_r + xi);
    Q['C'][2] = beta * pi_g;
    Q['C'][3] = alpha_1 * pi_t;

    Q['G'][0] = alpha_2 * pi_a;
    Q['G'][1] = beta * pi_c;
    Q['G'][2] = -(alpha_2 * pi_a + beta * pi_y + xi);
    Q['G'][3] = beta * pi_t;

    Q['T'][0] = beta * pi_a;
    Q['T'][1] = alpha_1 * pi_c;
    Q['T'][2] = beta * pi_g;
    Q['T'][3] = -(alpha_1 * pi_c + beta * pi_r + xi);

    return Q;
}


//' Initialize a MevoSampler object.
//'
//' @param Q An `unordered_map` of substitution rates for each nucleotide.
//' @param pis Vector of nucleotide equilibrium frequencies for
//'     "A", "C", "G", and "T", respectively.
//' @param xi Overall rate of indels.
//' @param psi Proportion of insertions to deletions.
//' @param rel_insertion_rates Relative insertion rates.
//' @param rel_deletion_rates Relative deletion rates.
//'
//' @noRd
//'
MevoSampler::MevoSampler(const std::unordered_map<char, std::vector<double>>& Q,
                         const double& xi,
                         const double& psi,
                         const std::vector<double>& pis,
                         arma::vec rel_insertion_rates,
                         arma::vec rel_deletion_rates) {

    uint n_events = 4 + rel_insertion_rates.n_elem + rel_deletion_rates.n_elem;

    // make relative rates sum to 1:
    rel_insertion_rates /= arma::accu(rel_insertion_rates);
    rel_deletion_rates /= arma::accu(rel_deletion_rates);
    // Now make them sum to the overall insertion/deletion rate:
    double xi_i = xi / (1 + 1/psi);  // overall insertion rate
    double xi_d = xi / (1 + psi);    // overall deletion rate
    rel_insertion_rates *= xi_i;
    rel_deletion_rates *= xi_d;

    /*
     (1) Combine substitution, insertion, and deletion rates into a single vector
     (2) Create TableSampler for each nucleotide
     (3) Fill the `rates` field with mutation rates for each nucleotide
     */
    double total_rate = 0;
    for (char c : mevo::bases) {
        std::vector<double> qc = Q.at(c);
        // Get the rate of change for this nucleotide
        double qi = -1.0 * qc[mevo::base_inds[c]];
        total_rate += qi;
        rates.q[c] = qi;
        /*
         Now that cell needs to be manually converted to zero so it's not sampled.
         (Remember, we want the probability of each, *given that a mutation occurs*.
          Mutating into itself doesn't count.)
         */
        qc[mevo::base_inds[c]] = 0;
        // Add insertions, then deletions
        qc.reserve(n_events);
        for (uint i = 0; i < rel_insertion_rates.n_elem; i++) {
            qc.push_back(rel_insertion_rates(i));
        }
        for (uint i = 0; i < rel_deletion_rates.n_elem; i++) {
            qc.push_back(rel_deletion_rates(i));
        }
        // Divide all by qi to make them probabilities
        for (uint i = 0; i < n_events; i++) qc[i] /= qi;
        // Fill TableSampler
        muts.sampler[c] = TableSampler(qc);
    }
    for (char c : mevo::bases) rates.w[c] = rates.q[c] / total_rate;

    // Now filling in event_lengths field of MutationSampler
    muts.event_lengths = std::vector<sint>(n_events, 0);
    for (uint i = 0; i < rel_insertion_rates.n_elem; i++) {
        muts.event_lengths[i + 4] = static_cast<sint>(i+1);
    }
    for (uint i = 0; i < rel_deletion_rates.n_elem; i++) {
        sint ds = static_cast<sint>(i + 1);
        ds *= -1;
        muts.event_lengths[i + 4 + rel_insertion_rates.n_elem] = ds;
    }

    // Now fill in the insertion-sequence sampler:
    nts = TableStringSampler<std::string>(table_sampler::bases, pis);

    return;
}




