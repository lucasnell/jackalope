
#include <RcppArmadillo.h>
#include <cmath>  // pow, log, exp
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class



#include "gemino_types.h"
#include "molecular_evolution.h"
#include "sequence_classes.h"  // Var* and Ref* classes
#include "pcg.h"  // pcg seeding
#include "table_sampler.h"  // table method of sampling

using namespace Rcpp;






//' Q matrix for rates for a given nucleotide using the TN93 substitution model.
//'
//' @noRd
//'
arma::mat TN93_rate_matrix(
        const double& pi_t, const double& pi_c,
        const double& pi_a, const double& pi_g,
        const double& alpha_1, const double& alpha_2, const double& beta,
        const double& xi) {

    arma::vec pis = {pi_t, pi_c, pi_a, pi_g};

    arma::mat Q(4,4);
    Q.fill(beta);
    Q.submat(arma::span(0,1), arma::span(0,1)).fill(alpha_1);
    Q.submat(arma::span(2,3), arma::span(2,3)).fill(alpha_2);
    for (uint i = 0; i < 4; i++) Q.col(i) *= pis(i);

    // Filling in diagonals
    Q.diag().fill(0.0);  // reset to zero so summing by row works
    arma::vec rowsums = arma::sum(Q, 1);
    rowsums += xi;
    rowsums *= -1;
    Q.diag() = rowsums;

    return Q;
}



//' Q matrix for rates for a given nucleotide using the JC69 substitution model.
//'
//' JC69 is a special case of TN93.
//'
//' @noRd
//'
arma::mat JC69_rate_matrix(
        const double& lambda, const double& xi) {

    arma::mat Q = TN93_rate_matrix(1, 1, 1, 1, lambda, lambda, lambda, xi);

    return Q;
}


//' Q matrix for rates for a given nucleotide using the K80 substitution model.
//'
//' K80 is a special case of TN93.
//'
//' @noRd
//'
arma::mat K80_rate_matrix(
        const double& alpha, const double& beta,
        const double& xi) {

    arma::mat Q = TN93_rate_matrix(1, 1, 1, 1, alpha, alpha, beta, xi);

    return Q;
}


//' Q matrix for rates for a given nucleotide using the F81 substitution model.
//'
//' F81 is a special case of TN93.
//'
//' @noRd
//'
arma::mat F81_rate_matrix(
        const double& pi_t, const double& pi_c,
        const double& pi_a, const double& pi_g,
        const double& xi) {

    arma::mat Q = TN93_rate_matrix(pi_a, pi_c, pi_g, pi_t, 1, 1, 1, xi);

    return Q;
}


//' Q matrix for rates for a given nucleotide using the HKY85 substitution model.
//'
//' HKY85 is a special case of TN93.
//'
//' @noRd
//'
arma::mat HKY85_rate_matrix(
        const double& pi_t, const double& pi_c,
        const double& pi_a, const double& pi_g,
        const double& alpha, const double& beta,
        const double& xi) {

    arma::mat Q = TN93_rate_matrix(pi_a, pi_c, pi_g, pi_t, alpha, alpha, beta, xi);


    return Q;
}


//' Q matrix for rates for a given nucleotide using the F84 substitution model.
//'
//' F84 is a special case of TN93.
//'
//' @noRd
//'
arma::mat F84_rate_matrix(
        const double& pi_t, const double& pi_c,
        const double& pi_a, const double& pi_g,
        const double& beta, const double& kappa,
        const double& xi) {

    double pi_y = pi_t + pi_c;
    double pi_r = pi_a + pi_g;

    double alpha_1 = 1 + kappa / pi_y;
    double alpha_2 = 1 + kappa / pi_r;

    arma::mat Q = TN93_rate_matrix(pi_a, pi_c, pi_g, pi_t, alpha_1, alpha_2, beta, xi);

    return Q;
}




//' Q matrix for rates for a given nucleotide using the GTR substitution model.
//'
//' @noRd
//'
arma::mat GTR_rate_matrix(
        const double& pi_t, const double& pi_c,
        const double& pi_a, const double& pi_g,
        const double& a, const double& b, const double& c,
        const double& d, const double& e, const double& f,
        const double& xi) {

    arma::vec pis = {pi_t, pi_c, pi_a, pi_g};

    arma::mat Q(4, 4, arma::fill::zeros);

    // Filling in non-diagonals
    arma::vec letters = {a, b, c, d, e, f};
    uint k = 0;
    for (uint i = 0; i < 3; i++) {
        for (uint j = i+1; j < 4; j++) {
            Q(i,j) = letters(k);
            Q(j,i) = letters(k);
            k++;
        }
    }
    for (uint i = 0; i < 4; i++) Q.col(i) *= pis(i);

    // Filling in diagonals
    arma::vec rowsums = arma::sum(Q, 1);
    rowsums += xi;
    rowsums *= -1;
    Q.diag() = rowsums;

    return Q;

}


//' Q matrix for rates for a given nucleotide using the UNREST substitution model.
//'
//' @param Q Matrix of rates for "T", "C", "A", and "G", respectively.
//'     Diagonal values are ignored.
//' @param xi Overall rate of indels.
//'
//' @noRd
//'
arma::mat UNREST_rate_matrix(
        arma::mat Q, const double& xi) {

    // reset to zero so summing by row works
    Q.diag().fill(0.0);
    // Filling in diagonals
    arma::vec rowsums = arma::sum(Q, 1);
    rowsums += xi;
    rowsums *= -1;
    Q.diag() = rowsums;

    return Q;

}





//' Initialize a MevoSampler object.
//'
//' @param Q A matrix of substitution rates for each nucleotide.
//' @param pis Vector of nucleotide equilibrium frequencies for
//'     "T", "C", "A", and "G", respectively.
//' @param xi Overall rate of indels.
//' @param psi Proportion of insertions to deletions.
//' @param rel_insertion_rates Relative insertion rates.
//' @param rel_deletion_rates Relative deletion rates.
//'
//' @noRd
//'
MevoSampler::MevoSampler(const arma::mat& Q,
                         const double& xi,
                         const double& psi,
                         const std::vector<double>& pis,
                         arma::vec rel_insertion_rates,
                         arma::vec rel_deletion_rates)
    : muts(), nts(), rates() {

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
    for (uint i = 0; i < 4; i++) {

        char c = mevo::bases[i];

        std::vector<double> qc = arma::conv_to<std::vector<double>>::from(Q.col(i));
        // Get the rate of change for this nucleotide
        double qi = -1.0 * qc[i];
        rates.rates[c] = qi;
        /*
         Now cell `i` in `qc` needs to be manually converted to zero so it's not sampled.
         (Remember, we want the probability of each, *given that a mutation occurs*.
          Mutating into itself doesn't count.)
         */
        qc[i] = 0;
        // Add insertions, then deletions
        qc.reserve(n_events);
        for (uint j = 0; j < rel_insertion_rates.n_elem; j++) {
            qc.push_back(rel_insertion_rates(j));
        }
        for (uint j = 0; j < rel_deletion_rates.n_elem; j++) {
            qc.push_back(rel_deletion_rates(j));
        }
        // Divide all by qi to make them probabilities
        for (uint j = 0; j < n_events; j++) qc[j] /= qi;
        // Fill TableSampler
        muts.sampler[i] = TableSampler(qc);
    }

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
    nts = TableStringSampler<std::string>(mevo::bases, pis);

    return;
}







/*
 At a time when an event occurs, this samples nucleotides based on their rates
and returns a random location where the event will occur

Reservoir sampling used from...
Efraimidis, P. S., and P. G. Spirakis. 2006. Weighted random sampling with a
reservoir. Information Processing Letters 97:181–185.
*/
uint event_location(const std::string& S,
                    const uint& chunk_size,
                    const MevoSampler& ms,
                    pcg32& eng) {

    // if (S.size() == 0) stop("Empty string sent to event_location.");
    // if (S.size() == 1) return 0;

    // Where should we start? Choose random location. (or 0 if chunk_size >= S.size())
    uint start;
    if (chunk_size < S.size()) {
        start = (static_cast<double>(eng()) / pcg::max) * (S.size() - chunk_size);
    } else start = 0;


    double r, key, X, w, t;
    // uint N = S.size();
    uint N = start + chunk_size;


    // Create a NucleoKeys object to store position and key
    r = runif_01(eng);
    key = std::pow(r, 1 / ms.rate(S[start]));
    double largest_key = key; // largest key (the one we're going to keep)
    // uint largest_pos = 0;     // position where largest key was found
    uint largest_pos = start;     // position where largest key was found

    // uint c = 0;
    uint c = start;
    while (c < (N-1)) {
        r = (static_cast<double>(eng()) + 1.0) / (pcg::max + 2.0);
        X = std::log(r) / std::log(largest_key);
        uint i = c + 1;
        double wt_sum0 = ms.rate(S[c]);
        double wt_sum1 = ms.rate(S[c]) + ms.rate(S[i]);
        while (X > wt_sum1 && i < (N-1)) {
            i++;
            wt_sum0 += ms.rate(S[(i-1)]);
            wt_sum1 += ms.rate(S[i]);
        }
        if (X > wt_sum1) break;
        if (wt_sum0 >= X) continue;

        largest_pos = i;

        w = ms.rate(S[i]);
        t = std::pow(largest_key, w);
        // r = runif_ab(eng, t, 1.0);
        r = t + ((static_cast<double>(eng()) + 1.0) / (pcg::max + 2.0)) * (1.0 - t);
        key = std::pow(r, 1 / w);
        largest_key = key;

        c = i;
    }

    return largest_pos;
}











//[[Rcpp::export]]
std::vector<uint> test_sampling(const std::string& seq, const uint& N,
                                const double& pi_t, const double& pi_c,
                                const double& pi_a, const double& pi_g,
                                const double& alpha_1, const double& alpha_2,
                                const double& beta,
                                const double& xi, const double& psi,
                                const arma::vec& rel_insertion_rates,
                                const arma::vec& rel_deletion_rates,
                                const uint& chunk_size) {

    arma::mat Q = TN93_rate_matrix(pi_t, pi_c, pi_a, pi_g, alpha_1, alpha_2, beta, xi);

    std::vector<double> pis = {pi_t, pi_c, pi_a, pi_g};

    MevoSampler ms(Q, xi, psi, pis, rel_insertion_rates, rel_deletion_rates);

    pcg32 eng = seeded_pcg();

    std::vector<uint> out(N);

    for (uint i = 0; i < N; i++) {
        Rcpp::checkUserInterrupt();
        out[i] = event_location(seq, chunk_size, ms, eng);
    }

    return out;
}



