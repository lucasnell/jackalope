#ifndef __JACKALOPE_MUTATOR_H
#define __JACKALOPE_MUTATOR_H


/*
 Combining samplers for substitutions and indels into a mutation sampler for
 evolving chromosomes along trees.
 */

#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <progress.hpp>  // for the progress bar
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // gamma_distribution


#include "jackalope_types.h" // integer types
#include "hap_classes.h"  // Hap* classes
#include "pcg.h"  // pcg seeding
#include "alias_sampler.h"  // alias method of sampling
#include "mutator_subs.h"   // SubMutator
#include "mutator_indels.h" // IndelMutator
#include "io.h"  // FileUncomp
#include "util.h"  // str_stop



using namespace Rcpp;





/*
 TreeMutator combines objects for adding substitutions and indels to chromosomes
 when evolving along trees.
 */
struct TreeMutator {

    // For adding substitutions:
    SubMutator subs;
    // For adding indels:
    IndelMutator indels;


    TreeMutator() {}

    TreeMutator(const std::vector<arma::mat>& Q_,
                    const std::vector<arma::mat>& U_,
                    const std::vector<arma::mat>& Ui_,
                    const std::vector<arma::vec>& L_,
                    const double& invariant_,
                    const arma::vec& insertion_rates,
                    const arma::vec& deletion_rates,
                    const double& epsilon,
                    const std::vector<double>& pi_tcag)
        : subs(Q_, U_, Ui_, L_, invariant_),
          indels(insertion_rates, deletion_rates, epsilon, pi_tcag) {}

    TreeMutator(const TreeMutator& other)
        : subs(other.subs), indels(other.indels) {}

    TreeMutator& operator=(const TreeMutator& other) {
        subs = other.subs;
        indels = other.indels;
        return *this;
    }

    /*
     Add mutations for a branch within a range.
     It also updates `end` for indels that occur in the range.
     (`end == begin` when chromosome region is of size zero bc `end` is non-inclusive)
     */
    int mutate(const double& b_len,
               HapChrom& hap_chrom,
               pcg64& eng,
               Progress& prog_bar,
               const uint64& begin,
               uint64& end,
               std::deque<uint8>& rate_inds);

    int new_rates(const uint64& begin,
                  const uint64& end,
                  std::deque<uint8>& rate_inds,
                  pcg64& eng,
                  Progress& prog_bar) {
        int status = subs.new_rates(begin, end, rate_inds, eng, prog_bar);
        return status;
    }

};












#endif
