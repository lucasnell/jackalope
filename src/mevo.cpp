
#include <RcppArmadillo.h>
#include <cmath>  // pow, log, exp
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <progress.hpp>  // for the progress bar



#include "gemino_types.hpp"
#include "mevo.hpp"
#include "seq_classes_var.hpp"  // Var* classes
#include "pcg.hpp"  // pcg seeding
#include "table_sampler.hpp"  // table method of sampling
#include "weighted_reservoir.hpp"  // weighted reservoir sampling
#include "mevo_gammas.hpp"  // SequenceGammas class

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
