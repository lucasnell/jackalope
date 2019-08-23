#ifndef __JACKALOPE_MUTATOR_INDELS_H
#define __JACKALOPE_MUTATOR_INDELS_H


#include "jackalope_types.h" // integer types and debugging preprocessor directives

/*
 This defines classes for adding insertions and deletions.
 */


#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // poisson_distribution


#include "var_classes.h"  // Var* classes
#include "pcg.h"  // pcg seeding
#include "alias_sampler.h"  // alias method of sampling
#include "util.h"  // str_stop
#include "mutator_subs.h"  // SubMutator



using namespace Rcpp;






class IndelMutator {

public:

    // Rates for each type of indel of each size, per bp per unit time:
    arma::vec rates;
    // Change in chromosome size associated with each indel (same length as `rates`):
    arma::vec changes;
    // Error control parameter (0 < `eps` << 1)
    double eps;
    // VarChrom object pointer to be manipulated:
    VarChrom* var_chrom;
    // For creating insertion sequences:
    AliasStringSampler<std::string> insert;

    IndelMutator() : var_chrom(nullptr) {}
    IndelMutator(const arma::vec& insertion_rates,
                 const arma::vec& deletion_rates,
                 const double& epsilon,
                 const std::vector<double>& pi_tcag)
        : rates(insertion_rates.n_elem + deletion_rates.n_elem),
          changes(insertion_rates.n_elem + deletion_rates.n_elem),
          eps(epsilon),
          insert("TCAG", pi_tcag),
          tau(0),
          rates_tau(insertion_rates.n_elem + deletion_rates.n_elem),
          n_events(insertion_rates.n_elem + deletion_rates.n_elem) {

        uint32 n = insertion_rates.n_elem;
        uint32 m = deletion_rates.n_elem;

        for (uint32 i = 0; i < n; i++) {
            rates(i) = insertion_rates(i);
            changes(i) = i+1;
        }
        for (uint32 i = 0; i < m; i++) {
            rates(i+n) = deletion_rates(i);
            changes(i+n) = -1.0 * static_cast<double>(i+1);
        }

    }

    IndelMutator(const IndelMutator& other)
        : rates(other.rates), changes(other.changes), eps(other.eps),
          var_chrom(other.var_chrom), insert(other.insert), tau(other.tau),
          rates_tau(other.rates_tau), n_events(other.n_events) {}

    IndelMutator& operator=(const IndelMutator& other) {
        rates = other.rates;
        changes = other.changes;
        eps = other.eps;
        var_chrom = other.var_chrom;
        insert = other.insert;
        tau = other.tau;
        rates_tau = other.rates_tau;
        n_events = other.n_events;
        return *this;
    }



    void new_chrom(VarChrom& var_chrom_) {
        var_chrom = &var_chrom_;
        return;
    }


    void add_indels(double b_len,
                    const uint64& begin,
                    const uint64& end,
                    SubMutator& subs,
                    pcg64& eng);


private:


    void calc_tau(double& b_len);

    // For generating # events per time period:
    std::poisson_distribution<uint32> distr = std::poisson_distribution<uint32>(1);
    // For "tau-leaping", `tau` is the time by which branch length can be split:
    double tau;
    // For storing rate over the entire chromosome over `tau` time units
    arma::vec rates_tau;
    // For storing # events that occur over `tau` time units, for each indel type + size
    std::vector<uint32> n_events;

};










#endif
