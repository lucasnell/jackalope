#ifndef __JACKALOPE_MUTATOR_INDELS_H
#define __JACKALOPE_MUTATOR_INDELS_H



/*
 This defines classes for adding insertions and deletions.
 */

#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <progress.hpp>  // for the progress bar
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // poisson_distribution


#include "jackalope_types.h" // integer types
#include "hap_classes.h"  // Hap* classes
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
    /*
     Error control parameter (0 <= `eps` << 1), where `eps==0` reverts to exact
     Doobs--Gillespie simulation
     */
    double eps;
    // For creating insertion sequences:
    AliasStringSampler<std::string> insert;
    // For total rate of all events if doing exact simulations
    double total_rate;
    // For sampling which event occurred if doing exact simulations
    AliasSampler event_sampler;

    IndelMutator() {}
    IndelMutator(const arma::vec& insertion_rates,
                 const arma::vec& deletion_rates,
                 const double& epsilon,
                 const std::vector<double>& pi_tcag)
        : rates(insertion_rates.n_elem + deletion_rates.n_elem),
          changes(insertion_rates.n_elem + deletion_rates.n_elem),
          eps(epsilon),
          insert("TCAG", pi_tcag),
          total_rate(0),
          event_sampler(),
          tau(0),
          rates_tau(insertion_rates.n_elem + deletion_rates.n_elem),
          n_events(insertion_rates.n_elem + deletion_rates.n_elem) {

        uint32 n = insertion_rates.n_elem;
        uint32 m = deletion_rates.n_elem;

        for (uint32 i = 0; i < n; i++) {
            rates(i) = insertion_rates(i);
            changes(i) = i+1;
            total_rate += insertion_rates(i);
        }
        for (uint32 i = 0; i < m; i++) {
            rates(i+n) = deletion_rates(i);
            changes(i+n) = -1.0 * static_cast<double>(i+1);
            total_rate += deletion_rates(i);
        }

        event_sampler = AliasSampler(rates.t());

    }

    IndelMutator(const IndelMutator& other)
        : rates(other.rates), changes(other.changes), eps(other.eps),
          insert(other.insert), total_rate(other.total_rate),
          event_sampler(other.event_sampler), tau(other.tau),
          rates_tau(other.rates_tau), n_events(other.n_events) {}

    IndelMutator& operator=(const IndelMutator& other) {
        rates = other.rates;
        changes = other.changes;
        eps = other.eps;
        insert = other.insert;
        total_rate = other.total_rate;
        event_sampler = other.event_sampler;
        tau = other.tau;
        rates_tau = other.rates_tau;
        n_events = other.n_events;
        return *this;
    }


    // Add indels, adjust `end` (`end == begin` when chromosome region is of size zero)
    int add_indels(double b_len,
                   const uint64& begin,
                   uint64& end,
                   std::deque<uint8>& rate_inds,
                   SubMutator& subs,
                   HapChrom& hap_chrom,
                   pcg64& eng,
                   Progress& prog_bar);


private:


    void calc_tau(double& b_len, HapChrom& hap_chrom);

    // For generating # events per time period:
    std::poisson_distribution<uint32> distr = std::poisson_distribution<uint32>(1);
    // For "tau-leaping", `tau` is the time by which branch length can be split:
    double tau;
    // For storing rate over the entire chromosome over `tau` time units
    arma::vec rates_tau;
    // For storing # events that occur over `tau` time units, for each indel type + size
    std::vector<uint32> n_events;
    // For generating jump lengths if doing exact simulations
    std::exponential_distribution<double> jump_distr =
        std::exponential_distribution<double>(1);


    // Add just one indel
    inline void one_indel__(const double& change,
                            std::string& insert_str,
                            const uint64& begin,
                            uint64& end,
                            std::deque<uint8>& rate_inds,
                            SubMutator& subs,
                            HapChrom& hap_chrom,
                            pcg64& eng);


    // Exact, rather than approximation to Doob--Gillespie algorithm
    inline void exact_sim(int& status,
                          double& b_len,
                          const uint64& begin,
                          uint64& end,
                          std::deque<uint8>& rate_inds,
                          SubMutator& subs,
                          HapChrom& hap_chrom,
                          pcg64& eng,
                          Progress& prog_bar);

};










#endif
