
#include <RcppArmadillo.h>
#include <cmath>  // pow, log, exp
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <progress.hpp>  // for the progress bar



#include "jackal_types.h"
#include "mevo.h"
#include "seq_classes_var.h"  // Var* classes
#include "pcg.h"  // pcg seeding
#include "table_sampler.h"  // table method of sampling
#include "weighted_reservoir.h"  // weighted reservoir sampling
#include "mevo_gammas.h"  // SequenceGammas class

using namespace Rcpp;





// Wrapper to make non-chunked version available from R

//[[Rcpp::export]]
SEXP make_mutation_sampler_base(const arma::mat& Q,
                                const std::vector<double>& pi_tcag,
                                const std::vector<double>& insertion_rates,
                                const std::vector<double>& deletion_rates) {

    XPtr<MutationSampler> out =
        make_mutation_sampler_base_<MutationSampler,LocationSampler>(
                Q, pi_tcag, insertion_rates, deletion_rates);

    return out;
}

// Same thing, but with chunks

//[[Rcpp::export]]
SEXP make_mutation_sampler_chunk_base(const arma::mat& Q,
                                      const std::vector<double>& pi_tcag,
                                      const std::vector<double>& insertion_rates,
                                      const std::vector<double>& deletion_rates,
                                      const uint32& chunk_size) {

    XPtr<ChunkMutationSampler> out =
        make_mutation_sampler_base_<ChunkMutationSampler,ChunkLocationSampler>(
                Q, pi_tcag, insertion_rates, deletion_rates);

    out->location.change_chunk(chunk_size);

    return out;
}





//' Used below to directly make a MutationTypeSampler
//'
//' @noRd
//'
MutationTypeSampler make_type_sampler(const arma::mat& Q,
                                      const std::vector<double>& pi_tcag,
                                      const std::vector<double>& insertion_rates,
                                      const std::vector<double>& deletion_rates) {

    std::vector<std::vector<double>> probs;
    std::vector<sint32> mut_lengths;
    std::vector<double> q_tcag;

    /*
    (1) Combine substitution, insertion, and deletion rates into a single vector
    (2) Fill the `q_tcag` vector with mutation rates for each nucleotide
    */
    fill_probs_q_tcag(probs, q_tcag, Q, pi_tcag, insertion_rates, deletion_rates);

    // Now filling in mut_lengths vector
    fill_mut_lengths(mut_lengths, insertion_rates, deletion_rates);

    // Type and insertion samplers:
    MutationTypeSampler type(probs, mut_lengths);

    return type;
}
