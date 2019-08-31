#ifndef __JACKALOPE_MUTATOR_H
#define __JACKALOPE_MUTATOR_H


/*
 Combining samplers for substitutions and indels into a mutation sampler for
 a single chromosome.
 */


#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // gamma_distribution


#include "jackalope_types.h" // integer types
#include "var_classes.h"  // Var* classes
#include "pcg.h"  // pcg seeding
#include "alias_sampler.h"  // alias method of sampling
#include "mutator_subs.h"   // SubMutator
#include "mutator_indels.h" // IndelMutator
#include "io.h"  // FileUncomp
#include "util.h"  // str_stop



using namespace Rcpp;





/*
 MutationSampler combines objects for sampling substitutions and indels.
 */
struct MutationSampler {

    // For adding substitutions:
    SubMutator subs;
    // For adding indels:
    IndelMutator indels;


    MutationSampler() {}

    MutationSampler(const std::vector<arma::mat>& Q_,
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

    MutationSampler(const MutationSampler& other)
        : subs(other.subs), indels(other.indels) {}

    MutationSampler& operator=(const MutationSampler& other) {
        subs = other.subs;
        indels = other.indels;
        return *this;
    }

    /*
     Add mutations for a branch within a range.
     It also updates `end` for indels that occur in the range.
     (`end == begin` when chromosome region is of size zero bc `end` is non-inclusive)
     */
    void mutate(const double& b_len, VarChrom& var_chrom, pcg64& eng,
                const uint64& begin, uint64& end)  {

#ifdef __JACKALOPE_DEBUG
        if (end < begin) stop("end < begin in MutationSampler.mutate");
        if (end == begin) stop("end == begin in MutationSampler.mutate");
#endif

        indels.add_indels(b_len, begin, end, subs, var_chrom, eng);

        subs.add_subs(b_len, begin, end, var_chrom, eng);

        return;
    }

};












#endif
