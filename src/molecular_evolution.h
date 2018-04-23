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
    const std::string bases = "ACGT";
    std::unordered_map<char, uint> base_inds = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
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
 */
class MutationRates {
public:

    std::unordered_map<char, double> q;  // rates for a given nucleotide
    std::unordered_map<char, double> w;  // same as above, but sums to 1

    MutationRates() {
        q = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};
        w = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};
    };
    MutationRates(const std::vector<double>& rates) {
        if (rates.size() != 4) {
            stop("Making a MutationRates object requires a vector of size 4");
        }
        q = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};
        w = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};
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
 Perhaps for a single mutation's info?
 */
struct MutationInfo {
    uint type;
    uint length;
};


/*
 For table-sampling events depending on which nucleotide you start with
 The event_types vector indicates which type of event a location in one of
     the unordered_map vectors is: 0:3 for substitution to A, C, G, and T,
     respectively, 4 for insertion, and 5 for deletion.
     It should be the same length as each vector for a given nucleotide in the
     `F` or `L` fields.
 The event_lengths vector should be the coincide with items in the event_types
     vector, but it tells how long each event is. For substitutions, this field
     is ignored, but a filler value should be provided so the event_lengths
     and event_types vectors are the same length.
 */
class EventSampler {
public:

    std::unordered_map<char, TableSampler> sampler;
    std::vector<sint> event_lengths;

    EventSampler() {
        sampler = {{'A', TableSampler()},
                   {'C', TableSampler()},
                   {'G', TableSampler()},
                   {'T', TableSampler()}};
    }
    EventSampler(const uint& N) : sampler(), event_lengths(N, 0) {
        sampler = {{'A', TableSampler()},
                   {'C', TableSampler()},
                   {'G', TableSampler()},
                   {'T', TableSampler()}};
    }
    // copy constructor
    EventSampler(const EventSampler& other)
        : sampler(other.sampler), event_lengths(other.event_lengths) {}
};



// MevoSampler combined objects for table-sampling event types and new nucleotides
class MevoSampler {
public:

    MevoSampler(const EventSampler& event_, const TableStringSampler<std::string>& nucleo_)
        : event_sampler(event_), nt_sampler(nucleo_) {};
    // MevoSampler(uint N) : event(N), nucleo() {};
    // MevoSampler() : event(), nucleo() {};
private:
    EventSampler event_sampler;
    TableStringSampler<std::string> nt_sampler;
};


#endif
