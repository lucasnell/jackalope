#ifndef __GEMINO_MEVO_H
#define __GEMINO_MEVO_H


#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // gamma_distribution
#include <unordered_map>  // unordered_map


#include "gemino_types.h" // integer types
#include "sequence_classes.h"  // Var* and Ref* classes
#include "pcg.h"  // pcg seeding
#include "table_sampler.h"  // table method of sampling



using namespace Rcpp;

namespace mevo {
    const std::string bases = "TCAG";
    std::unordered_map<char, uint> base_inds = {{'T', 0}, {'C', 1}, {'A', 2}, {'G', 3}};
}


/*
 One sequence's Gamma values
 THIS IS TOO SIMPLE: INDELS CHANGE RATES OF SEQUENCES OUTSIDE THE INITIAL RANGE
 ... perhaps unless you use old positions rather than new ones to find Gamma values
 */
struct SeqGammas {
    std::vector<double> gamma_vals;  // values from a Gamma distribution
    uint k;                          // number of bp per value

    SeqGammas(const uint& k_, const RefSequence& rs,
              pcg32& eng, const double& alpha) : gamma_vals(), k(k_) {
        // Now fill gammas vector
        uint n_gammas = static_cast<uint>(std::ceil(
            static_cast<double>(rs.size()) / static_cast<double>(k_)));
        gamma_vals = std::vector<double>(n_gammas);
        std::gamma_distribution<double> distr(alpha, alpha);
        for (uint i = 0; i < n_gammas; i++) {
            gamma_vals[i] = distr(eng);
        }
    }

    double get_gamma(const uint& pos) {
        uint ind = pos / k;
        if (ind >= gamma_vals.size()) stop("ind too high for this SeqGammas object");
        return gamma_vals[ind];
    }
};




/*
 Stores info on the overall mutation rates for each nucleotide.
 Ns are set to 0 bc we don't want to process these.
 */
class MutationRates {
public:

    std::unordered_map<char, double> q;  // rates for a given nucleotide
    std::unordered_map<char, double> w;  // same as above, but sums to 1

    MutationRates() {
        q = {{'T', 0}, {'C', 0}, {'A', 0}, {'G', 0}, {'N', 0}};
        w = {{'T', 0}, {'C', 0}, {'A', 0}, {'G', 0}, {'N', 0}};
    };
    MutationRates(const std::vector<double>& rates) {
        if (rates.size() != 4) {
            stop("Making a MutationRates object requires a vector of size 4");
        }
        q = {{'T', 0}, {'C', 0}, {'A', 0}, {'G', 0}, {'N', 0}};
        w = {{'T', 0}, {'C', 0}, {'A', 0}, {'G', 0}, {'N', 0}};
        for (uint i = 0; i < 4; i++) q[mevo::bases[i]] = rates[i];
        double rate_sum = std::accumulate(rates.begin(), rates.end(), 0.0);
        for (uint i = 0; i < 4; i++) w[mevo::bases[i]] = rates[i] / rate_sum;
    };

};




// struct to store info on a nucleotide's key value and position for weighted sampling
struct NucleoKeyPos {

    double key;
    uint pos;

    NucleoKeyPos() : key(0.0), pos(0) {};
    NucleoKeyPos(const double& _key, const uint& _pos) : key(_key), pos(_pos) {};

};


/*
 For a single mutation's info
 */
struct MutationInfo {
    char nucleo;
    sint length;
    // Initialize from an index and event-lengths vector
    MutationInfo (const uint& ind, const std::vector<sint>& event_lengths)
        : nucleo('\0'), length(0) {
        if (ind < 4) {
            nucleo = mevo::bases[ind];
        } else {
            length = event_lengths[ind];
        }
    }
};


/*
 For table-sampling events depending on which nucleotide you start with.
 The event_lengths vector tells how long each event is.
 This field is 0 for substitions, < 0 for deletions, and > 0 for insertions.
 */
class MutationSampler {
public:

    std::unordered_map<char, TableSampler> sampler;
    std::vector<sint> event_lengths;

    MutationSampler() {
        sampler = {{'A', TableSampler()},
                   {'C', TableSampler()},
                   {'G', TableSampler()},
                   {'T', TableSampler()}};
    }
    MutationSampler(const uint& N) : sampler(), event_lengths(N, 0) {
        sampler = {{'A', TableSampler()},
                   {'C', TableSampler()},
                   {'G', TableSampler()},
                   {'T', TableSampler()}};
    }
    // copy constructor
    MutationSampler(const MutationSampler& other)
        : sampler(other.sampler), event_lengths(other.event_lengths) {}

    /*
     Sample an event based on an input nucleotide.
     Will return error if `c` is not one of 'A', 'C', 'G', or 'T'
     So make sure no 'N's get input here!
     */
    MutationInfo sample(const char& c, pcg32& eng) const {
        uint ind = sampler.at(c).sample(eng);
        MutationInfo mi(ind, event_lengths);
        return mi;
    }
};



// MevoSampler combined objects for table-sampling event types and new nucleotides
class MevoSampler {
public:

    MevoSampler(const std::unordered_map<char, std::vector<double>>& Q,
                const double& xi, const double& psi, const std::vector<double>& pis,
                arma::vec rel_insertion_rates, arma::vec rel_deletion_rates);
    MevoSampler(const MutationSampler& event_,
                const TableStringSampler<std::string>& nucleo_,
                const MutationRates& rates_)
        : muts(event_), nts(nucleo_), rates(rates_) {};

private:
    // For sampling the type of mutation:
    MutationSampler muts;
    // For insertion sequences:
    TableStringSampler<std::string> nts;
    // For overall mutation rates by nucleotide:
    MutationRates rates;
};


#endif
