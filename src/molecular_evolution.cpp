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











//' Q matrix for rates for a given nucleotide using the TN93 substitution model.
//'
//' @noRd
//'
std::unordered_map<char, std::vector<double>> TN93_rate_matrix(
        const double& pi_t, const double& pi_c,
        const double& pi_a, const double& pi_g,
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

    Q['T'] = {-(alpha_1 * pi_c + beta * pi_r + xi),
              alpha_1 * pi_c,
              beta * pi_a,
              beta * pi_g};

    Q['C'] = {alpha_1 * pi_t,
              -(alpha_1 * pi_t + beta * pi_r + xi),
              beta * pi_a,
              beta * pi_g};

    Q['A'] = {beta * pi_t,
              beta * pi_c,
              -(alpha_2 * pi_g + beta * pi_y + xi),
              alpha_2 * pi_g};

    Q['G'] = {beta * pi_t,
              beta * pi_c,
              alpha_2 * pi_a,
              -(alpha_2 * pi_a + beta * pi_y + xi)};

    return Q;
}



//' Q matrix for rates for a given nucleotide using the JC69 substitution model.
//'
//' JC69 is a special case of TN93.
//'
//' @noRd
//'
std::unordered_map<char, std::vector<double>> JC69_rate_matrix(
        const double& lambda, const double& xi) {

    std::unordered_map<char, std::vector<double>> Q;
    Q = TN93_rate_matrix(1, 1, 1, 1, lambda, lambda, lambda, xi);

    return Q;
}


//' Q matrix for rates for a given nucleotide using the K80 substitution model.
//'
//' K80 is a special case of TN93.
//'
//' @noRd
//'
std::unordered_map<char, std::vector<double>> K80_rate_matrix(
        const double& alpha, const double& beta,
        const double& xi) {

    std::unordered_map<char, std::vector<double>> Q;
    Q = TN93_rate_matrix(1, 1, 1, 1, alpha, alpha, beta, xi);

    return Q;
}


//' Q matrix for rates for a given nucleotide using the F81 substitution model.
//'
//' F81 is a special case of TN93.
//'
//' @noRd
//'
std::unordered_map<char, std::vector<double>> F81_rate_matrix(
        const double& pi_t, const double& pi_c,
        const double& pi_a, const double& pi_g,
        const double& xi) {

    std::unordered_map<char, std::vector<double>> Q;
    Q = TN93_rate_matrix(pi_a, pi_c, pi_g, pi_t, 1, 1, 1, xi);

    return Q;
}


//' Q matrix for rates for a given nucleotide using the HKY85 substitution model.
//'
//' HKY85 is a special case of TN93.
//'
//' @noRd
//'
std::unordered_map<char, std::vector<double>> HKY85_rate_matrix(
        const double& pi_t, const double& pi_c,
        const double& pi_a, const double& pi_g,
        const double& alpha, const double& beta,
        const double& xi) {

    std::unordered_map<char, std::vector<double>> Q;
    Q = TN93_rate_matrix(pi_a, pi_c, pi_g, pi_t, alpha, alpha, beta, xi);


    return Q;
}


//' Q matrix for rates for a given nucleotide using the F84 substitution model.
//'
//' F84 is a special case of TN93.
//'
//' @noRd
//'
std::unordered_map<char, std::vector<double>> F84_rate_matrix(
        const double& pi_t, const double& pi_c,
        const double& pi_a, const double& pi_g,
        const double& beta, const double& kappa,
        const double& xi) {

    double pi_y = pi_t + pi_c;
    double pi_r = pi_a + pi_g;

    double alpha_1 = 1 + kappa / pi_y;
    double alpha_2 = 1 + kappa / pi_r;

    std::unordered_map<char, std::vector<double>> Q;
    Q = TN93_rate_matrix(pi_a, pi_c, pi_g, pi_t, alpha_1, alpha_2, beta, xi);

    return Q;
}




//' Q matrix for rates for a given nucleotide using the GTR substitution model.
//'
//' @noRd
//'
std::unordered_map<char, std::vector<double>> GTR_rate_matrix(
        const double& pi_t, const double& pi_c,
        const double& pi_a, const double& pi_g,
        const double& a, const double& b, const double& c,
        const double& d, const double& e, const double& f,
        const double& xi) {

    arma::mat Qmat(4, 4, arma::fill::zeros);

    // Filling in non-diagonals
    arma::vec letters = {a, b, c, d, e, f};
    uint k = 0;
    for (uint i = 0; i < 3; i++) {
        for (uint j = i+1; j < 4; j++) {
            Qmat(i,j) = letters(k);
            Qmat(j,i) = letters(k);
            k++;
        }
    }
    Qmat.col(0) *= pi_t;
    Qmat.col(1) *= pi_c;
    Qmat.col(2) *= pi_a;
    Qmat.col(3) *= pi_g;

    // Filling in diagonals
    arma::vec rowsums = arma::sum(Qmat, 1);
    rowsums += xi;
    rowsums *= -1;
    Qmat.diag() = rowsums;

    // Converting to unordered_map
    std::unordered_map<char, std::vector<double>> Q;
    Q = {{'T', arma::conv_to<std::vector<double>>::from(Qmat.row(0))},
         {'C', arma::conv_to<std::vector<double>>::from(Qmat.row(1))},
         {'A', arma::conv_to<std::vector<double>>::from(Qmat.row(2))},
         {'G', arma::conv_to<std::vector<double>>::from(Qmat.row(3))}};

    return Q;

}


//' Q matrix for rates for a given nucleotide using the UNREST substitution model.
//'
//' @param Qmat Matrix of rates for "T", "C", "A", and "G", respectively.
//'     Diagonals are ignored.
//' @param xi Overall rate of indels.
//'
//' @noRd
//'
std::unordered_map<char, std::vector<double>> UNREST_rate_matrix(
        arma::mat Qmat, const double& xi) {

    // Filling in diagonals
    arma::vec rowsums = arma::sum(Qmat, 1);
    rowsums += xi;
    rowsums *= -1;
    Qmat.diag() = rowsums;

    // Converting to unordered_map
    std::unordered_map<char, std::vector<double>> Q;
    Q = {{'T', arma::conv_to<std::vector<double>>::from(Qmat.row(0))},
    {'C', arma::conv_to<std::vector<double>>::from(Qmat.row(1))},
    {'A', arma::conv_to<std::vector<double>>::from(Qmat.row(2))},
    {'G', arma::conv_to<std::vector<double>>::from(Qmat.row(3))}};

    return Q;

}





//' Initialize a MevoSampler object.
//'
//' @param Q An `unordered_map` of substitution rates for each nucleotide.
//' @param pis Vector of nucleotide equilibrium frequencies for
//'     "T", "C", "A", and "G", respectively.
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







/*
 At a time when an event occurs, this samples nucleotides based on their rates
and returns a random location where the event will occur

Reservoir sampling used from...
Efraimidis, P. S., and P. G. Spirakis. 2006. Weighted random sampling with a
reservoir. Information Processing Letters 97:181–185.
*/
uint event_location(const std::string& S,
                    const MevoSampler& ms,
                    pcg32& eng) {

    if (S.size() == 0) stop("Empty string sent to event_location.");
    if (S.size() == 1) return 0;

    double r, key, X, w, t;
    uint N = S.size();

    // Create a NucleoKeys object to store position and key
    r = runif_01(eng);
    key = std::log(r) / ms.rate(S[0]);
    key = std::exp(key);
    NucleoKeyPos pq(key, 0);

    uint c = 0;
    while (c < (N-1)) {
        r = runif_01(eng);
        X = std::log(r) / std::log(pq.key);
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

        pq.pos = i;

        w = ms.rate(S[i]);
        t = std::pow(pq.key, w);
        r = runif_ab(eng, t, 1.0);
        key = std::pow(r, 1 / w);
        pq.key = key;

        c = i;
    }

    return pq.pos;
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
                                const uint& print_every = 1000) {

    std::unordered_map<char, std::vector<double>> Q;
    Q = TN93_rate_matrix(pi_t, pi_c, pi_a, pi_g, alpha_1, alpha_2, beta, xi);

    for (const char& c : mevo::bases) {
        for (const double& x : Q.at(c)) Rcout << x << ' ';
        Rcout << std::endl;
    }

    std::vector<double> pis = {pi_t, pi_c, pi_a, pi_g};

    MevoSampler ms(Q, xi, psi, pis, rel_insertion_rates, rel_deletion_rates);

    pcg32 eng = seeded_pcg();

    std::vector<uint> out(N);

    for (uint i = 0; i < N; i++) {
        Rcpp::checkUserInterrupt();
        out[i] = event_location(seq, ms, eng);
        if (i % print_every == 0) Rcout << i << std::endl;
    }

    return out;
}



