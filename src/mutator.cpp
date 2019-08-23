

#include "jackalope_types.h"

#include <RcppArmadillo.h>
#include <cmath>  // exp
#include <vector>  // vector class
#include <string>  // string class


#include "mutator_subs.h"   // SubMutator
#include "mutator_indels.h" // IndelMutator
#include "mutator.h"

using namespace Rcpp;






/*
 Add mutations for a branch within a range.
 It also updates `end` for indels that occur in the range.
 (`end == begin` when chromosome region is of size zero bc `end` is non-inclusive)
 */
void MutationSampler::mutate(const double& b_len,
                             pcg64& eng,
                             const uint64& begin,
                             uint64& end) {

#ifdef __JACKALOPE_DEBUG
    if (end < begin) stop("end < begin in MutationSampler.mutate");
    if (end == begin) stop("end == begin in MutationSampler.mutate");
#endif

    indels.add_indels(b_len, begin, end, subs, eng);

    subs.add_subs(b_len, begin, end, eng);

    return;
}





// Wrapper to make mutation sampler available from R

//[[Rcpp::export]]
SEXP make_mutation_sampler_base(const std::vector<arma::mat>& Q,
                                const std::vector<arma::mat>& U,
                                const std::vector<arma::mat>& Ui,
                                const std::vector<arma::vec>& L,
                                const double& invariant,
                                const arma::vec& insertion_rates,
                                const arma::vec& deletion_rates,
                                const double& epsilon,
                                const std::vector<double>& pi_tcag) {


    /*
     Now create and fill output pointer to base sampler:
     */
    XPtr<MutationSampler> mutator(new MutationSampler());

    mutator->subs = SubMutator(Q, U, Ui, L, invariant);
    mutator->indels = IndelMutator(insertion_rates, deletion_rates, epsilon, pi_tcag);

    return mutator;
}


