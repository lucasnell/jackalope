
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
#include "weighted_reservoir.h"  // weighted reservoir sampling
#include "mevo_gammas.h"  // SequenceGammas class
#include "mevo_rate_matrices.h"  // rate matrix functions

using namespace Rcpp;





//' Initialize a MutationSampler object.
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
MutationSampler::MutationSampler(const arma::mat& Q,
                                 const double& xi,
                                 const double& psi,
                                 const std::vector<double>& pis,
                                 arma::vec rel_insertion_rates,
                                 arma::vec rel_deletion_rates)
    : rates(), types(), nucleos() {

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
        rates[c] = qi;
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
        types.sampler[i] = TableSampler(qc);
    }

    // Now filling in event_lengths field of MutationTypeSampler
    types.event_lengths = std::vector<sint>(n_events, 0);
    for (uint i = 0; i < rel_insertion_rates.n_elem; i++) {
        types.event_lengths[i + 4] = static_cast<sint>(i+1);
    }
    for (uint i = 0; i < rel_deletion_rates.n_elem; i++) {
        sint ds = static_cast<sint>(i + 1);
        ds *= -1;
        types.event_lengths[i + 4 + rel_insertion_rates.n_elem] = ds;
    }

    // Now fill in the insertion-sequence sampler:
    nucleos = TableStringSampler<std::string>(mevo::bases, pis);

    return;
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

    MutationSampler ms(Q, xi, psi, pis, rel_insertion_rates, rel_deletion_rates);

    pcg32 eng = seeded_pcg();

    std::vector<uint> out(N);

    for (uint i = 0; i < N; i++) {
        Rcpp::checkUserInterrupt();
        out[i] = event_location(seq, chunk_size, ms.rates, eng);
    }

    return out;
}



